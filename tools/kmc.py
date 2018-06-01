'''
    Tool wrapper for the KMC kmer counter 
    ( http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about )
'''

import logging
import os
import os.path
import subprocess
import shutil
import shlex

import tools
import tools.samtools
import tools.picard
import util.file
import util.misc

import Bio.SeqIO

TOOL_NAME = 'kmc'
TOOL_VERSION = '3.1.0'

log = logging.getLogger(__name__)

MAX_COUNT=999999  # an int larger than any realistic kmer count or read length

class KmcTool(tools.Tool):
    '''Tool wrapper for KMC kmer counter'''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, executable='kmc',
                                                  verifycmd='kmc -h > /dev/null')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    @staticmethod
    def _kmc_db_name(kmc_db):
        """Return the kmer database path, given either the db path of the file name of the .kmc_pre or .kmc_suf file"""
        base, ext = os.path.splitext(kmc_db)
        return base if ext in ('.kmc_pre', '.kmc_suf') else kmc_db

    @staticmethod
    def _get_file_format_opt(fname):
        """Get the KMC command-line option to specify file format"""
        file_type = util.file.uncompressed_file_type(fname)
        if file_type in ('.fasta', '.fa'): return 'm'
        elif file_type in ('.fastq', '.fq'): return 'q'
        elif file_type == '.bam': return 'bam'
        else: raise RuntimeError('Unknown seq file format: {}'.format(file_type))

    def build_kmer_db(self, seq_files, kmer_size, min_occs, max_occs, counter_cap, kmc_db, mem_limit_gb=8, threads=None, kmc_opts=''):
        """Build a database of kmers occurring in the given sequence files.

        Inputs:
          seq_files: filename(s) of files from which to gather kmers.  Files can be fasta, fastq or bam; fasta and fastq
             may be compressed with gzip or bzip2.

        Params:
          kmer_size: kmer size
          min_occs: ignore kmers occurring fewer than this many times
          max_occs: ignore kmers occurring more than this many times
          counter_cap: when writing kmer counts to the database, cap the values at this number
          kmc_opts: any additional kmc flags

        Outputs:
          kmc_db: kmc database path, with or without the .kmc_pre/.kmc_suf extension.
          
        """
        seq_files = util.misc.make_seq(seq_files)
        kmc_db = self._kmc_db_name(kmc_db)
        threads = util.misc.sanitize_thread_count(threads)
        with util.file.tmp_dir(suffix='kmcdb') as t_dir, util.file.tempfname(suffix='kmcfiles') as seq_file_list:
            util.file.dump_file(seq_file_list, '\n'.join(seq_files))
            args = ['-v']

            input_fmt_opts = set(map(self._get_file_format_opt, seq_files))
            assert len(input_fmt_opts) == 1, "All input files must be of the same format"
            input_fmt_opt = list(input_fmt_opts)[0]

            args += '-f{} -k{} -ci{} -cx{} -cs{} -m{} -sm -t{} @{} {}'.format(input_fmt_opt, kmer_size, min_occs, max_occs, counter_cap, mem_limit_gb, threads,
                                                                              seq_file_list, kmc_db).split()
            if kmc_opts: args += shlex.split(kmc_opts)
            args += [t_dir]
            tool_cmd = [self.install_and_get_path()] + args
            log.info('Building KMC database with command: ' + ' '.join(tool_cmd))
            subprocess.check_call(tool_cmd)
            assert os.path.isfile(kmc_db+'.kmc_pre') and os.path.isfile(kmc_db+'.kmc_suf'), 'KMC database files not created: {}'.format(kmc_db)

    def execute(self, args, threads=None):
        """Run kmc_tools with the given args"""
        threads = util.misc.sanitize_thread_count(threads)
        tool_cmd = [self.install_and_get_path()+'_tools'] + ['-v', '-t{}'.format(threads)] + list(map(str, args))
        log.info('Running kmc_tools command: ' + ' '.join(tool_cmd))
        print('Running kmc_tools command: ' + ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)


    def dump_kmers(self, kmc_db, out_kmers, min_occs=1, max_occs=MAX_COUNT, threads=None):
        """Dump the kmers from the database to a text file"""
        self.execute('transform {} -ci{} -cx{} dump -s {}'.format(self._kmc_db_name(kmc_db), min_occs, max_occs, out_kmers).split(), threads=threads)
        assert os.path.isfile(out_kmers)


    def filter_reads(self, kmc_db, in_reads, out_reads, db_min_occs=1, db_max_occs=MAX_COUNT, reads_min_occs=None, reads_max_occs=None, threads=None):
        """Filter reads based on their kmer contents.

        Note that "occurrence of a kmer" means "occurrence of the kmer or its reverse complement".

        Inputs:
          kmc_db: the kmc kmer database
          in_reads: the reads to filter.  can be a .fasta or .fastq or .bam; fasta or fastq can be compressed with gzip or bzip2.
             If a .bam, a read pair is kept if either mate passes the filter.

        Outputs:
          out_reads: file to which filtered reads are written.  type is determined from extension, same types as above are supported.

        Params:
          db_min_occs: only consider database kmers with at least this count
          db_max_occs: only consider database kmers with at most this count
          reads_min_occs: only keep reads with at least this many occurrences of kmers from database.  If `as_perc` is True, interpreted as percent
             of read length
          reads_max_occs: only keep reads with no more than this many occurrence of kmers from the database
        """

        if reads_min_occs is None: reads_min_occs = 0
        if reads_max_occs is None: reads_max_occs = MAX_COUNT

        print('IN FILER_READS')

        threads = util.misc.sanitize_thread_count(threads)
        in_reads_type = util.file.uncompressed_file_type(in_reads)
        _in_reads = in_reads
        _out_reads = out_reads
        with util.file.tmp_dir(suffix='kmcfilt') as t_dir:
            if in_reads_type in ('.fa', '.fasta'):
                # kmc_tools filter currently requires fasta files to be in fasta-2line format
                # https://github.com/refresh-bio/KMC/issues/57
                _in_reads = os.path.join(t_dir, 'in_reads.fasta')
                with util.file.open_or_gzopen(in_reads, 'rt'):
                    Bio.SeqIO.convert(in_reads, 'fasta', _in_reads, 'fasta-2line')
            if in_reads_type == '.bam':
                # kmc_tools filter currently does not support .bam files
                # https://github.com/refresh-bio/KMC/issues/66
                _in_reads = os.path.join(t_dir, 'in_reads.fasta')
                tools.samtools.SamtoolsTool().bam2fa(in_reads, _in_reads, add_mate_num=True)
                _out_reads = os.path.join(t_dir, 'out_reads.fasta')

            in_reads_fmt = 'q' if in_reads_type in ('.fq', '.fastq') else 'a'

            self.execute('filter {} -ci{} -cx{} {} -f{} -ci{} -cx{} {}'.format(self._kmc_db_name(kmc_db), db_min_occs, db_max_occs, _in_reads, in_reads_fmt, 
                                                                               reads_min_occs, reads_max_occs, _out_reads).split(), threads=threads)

            if in_reads_type == '.bam':
                assert out_reads.endswith('.bam')
                passed_read_names = os.path.join(t_dir, 'passed_read_names.txt')
                self._get_fasta_read_names(_out_reads, passed_read_names)
                tools.picard.FilterSamReadsTool().execute(inBam=in_reads, exclude=False, readList=passed_read_names, outBam=out_reads)
        # end: with util.file.tmp_dir(suffix='kmcfilt') as t_dir
    # end: def filter_reads(self, kmc_db, in_reads, out_reads, db_min_occs=1, db_max_occs=MAX_COUNT, reads_min_occs=None, reads_max_occs=None, threads=None)

    @staticmethod
    def _get_fasta_read_names(in_fasta, out_read_names):
        """Save the read names of reads in a .fasta file to a text file"""
        with open(in_fasta) as in_fasta_f, open(out_read_names, 'wt') as out_read_names_f:
            for line in in_fasta_f:
                if line.startswith('>'):
                    assert line.endswith('/1\n') or line.endswith('/2\n')
                    out_read_names_f.write(line[1:-3]+'\n')

# end: class KmcTool(tools.Tool)

