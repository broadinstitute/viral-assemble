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
import util.file
import util.misc

import Bio.SeqIO

TOOL_NAME = 'kmc'
TOOL_VERSION = '3.1.0'

log = logging.getLogger(__name__)

class KmcTool(tools.Tool):
    '''Tool wrapper for KMC kmer counter'''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, executable='kmc',
                                                  verifycmd='kmc -h > /dev/null')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def kmc_db_name(self, kmc_db):
        """Return the kmer database path, given either the db path of the file name of the .kmc_pre or .kmc_suf file"""
        base, ext = os.path.splitext(kmc_db)
        return base if ext in ('.kmc_pre', '.kmc_suf') else kmc_db
        

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
        kmc_db = self.kmc_db_name(kmc_db)
        threads = util.misc.sanitize_thread_count(threads)
        with util.file.tmp_dir(suffix='kmcdb') as t_dir, util.file.tempfname(suffix='kmcfiles') as seq_file_list:
            util.file.dump_file(seq_file_list, '\n'.join(seq_files))
            args = ['-v']

            def get_input_fmt_opt(fname):
                """Get the KMC command-line option to specify input file format"""
                fname_base, fname_ext = os.path.splitext(fname)
                if fname_ext in ('.gz', '.bz2'):
                    fname_ext = os.path.splitext(fnamebase)[1]
                if fname_ext in ('.fasta', '.fa'): return 'm'
                elif fname_ext in ('.fastq', '.fq'): return 'q'
                elif fname_ext == '.bam': return 'bam'
                else: raise RuntimeError('Unknown seq file format: {}'.format(fname_ext))

            input_fmt_opts = set(map(get_input_fmt_opt, seq_files))
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

    def execute(self, args):
        tool_cmd = [self.install_and_get_path()+'_tools'] + list(map(str, args))
        log.info('Running kmc_tools command: ' + ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def dump_kmers(self, kmc_db, out_kmers, min_occs=1, max_occs=999999, threads=None):
        """Dump the kmers from the database to a text file"""
        threads = util.misc.sanitize_thread_count(threads)
        self.execute('-t{} transform {} -ci{} -cx{} dump -s {}'.format(threads, self.kmc_db_name(kmc_db), min_occs, max_occs, out_kmers).split())
        assert os.path.isfile(out_kmers)

# end: class KmcTool(tools.Tool)

