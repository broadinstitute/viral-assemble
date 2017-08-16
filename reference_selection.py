#!/usr/bin/env python
'''Tools for selecting, for each sample, the reference closest to that sample.  "refsel" refers to "reference
selection" throughout.
'''

__author__ = "ilya@broadinstitute.org"
__commands__ = []

# built-ins
import argparse
import collections
import logging
import os
import os.path
import shutil

# intra-module
import util.cmd
import util.file
import util.misc
import tools.clark

# third-party
import Bio.SeqIO
import Bio.Seq
import Bio.SeqRecord
import Bio.Alphabet

_log = logging.getLogger(__name__)   # pylint: disable=C0103

_DB_FORMAT_VERSION = "1.1"

def _db_format_string():
    '''Return a string identifying the format of the refsel database.'''
    return '|'.join(('RefSelDB_'+ _DB_FORMAT_VERSION, tools.clark.TOOL_NAME, tools.clark.TOOL_VERSION))

class _RefSelPaths(object):  # pylint: disable=R0902
    '''File and directory paths used by the refsel module.'''

    def __init__(self, refsel_dir):
        '''Construct all the paths used by refsel code, relative to `refsel-dir`'''

        from os.path import join

        self.refsel_dir = refsel_dir

        self.db_format = join(refsel_dir, 'db_format.dat')
        self.refs_fasta = join(refsel_dir, 'refs.fasta')
        self.n_segments_per_genome = join(refsel_dir, 'n_segments_per_genome.dat')
        self.kmer_sizes = join(refsel_dir, 'kmer_sizes.dat')

        self.refs_segments_dir = join(refsel_dir, 'refs_segments')

        self.clark_dir = join(refsel_dir, 'CLARK')
        self.clark_variant = join(self.clark_dir, 'clark_variant')
        self.clark_targets_addresses = join(self.clark_dir, 'targets_addresses.txt')

    # end: def __init__(self, refsel_dir)

    def ref_segment_fasta(self, ref_id, segment_id):
        '''Filename of fasta file containing one segment of the given reference.

        Args:
            ref_id: our internal ID of the given reference (we use the FASTA id of the first segment of a reference
                as reference ID)
            segment_id: sequence ID for the segment
        '''
        return os.path.join(self.refs_segments_dir,
                            util.file.string_to_file_name(ref_id)+'_'+util.file.string_to_file_name(segment_id)+'.fa')

    def clark_db_dir(self, kmer_size):
        '''Path to CLARK database directory for given kmer size'''
        return os.path.join(self.clark_dir, 'k'+str(kmer_size))

# end: class _RefSelPaths(object)

# =========================
# ***  build_refsel_db  ***
# =========================

def build_refsel_db(refs_segments, refs_defs, n_segments_per_genome, refsel_dir, # pylint: disable=R0914
                    kmer_sizes=None, clark_variant='light'):
    '''Construct the refsel database.  The database is written to the directory `refsel_dir` .

    Args:
        refs_segments: fasta file giving the nucleotide sequences of the segments of the reference genomes.
        refs_defs: text file defining the segments comprising each reference, and 
        n_segments_per_genome: number of segments in the genome.  The number of fasta records in `refs_fasta`
            must be a multiple of this number.
        refsel_dir: directory under which the refsel database will be stored
        kmer_sizes: kmer size(s) for which to build the database.  If not given, defaults for (31,) for CLARK and
            (27,) for CLARK-l (the low-memory, less-precise 'light' variant of CLARK).  Only the first kmer size is 
            currently used for reference selection, but this may change in the future.
        clark_variant: which variant of the CLARK metagenomics tool to use.
            Options are 'standard' and 'light' (default).  The light variant uses less memory
            but is less precise.
    '''

    assert clark_variant in ('standard', 'light')
    _sanity_check_ref_genomes(refs_fasta, n_segments_per_genome)

    from util.file import mkdir_p, dump_file

    paths = _RefSelPaths(refsel_dir)
    del refsel_dir  # make sure all paths are accessed via the `paths` object

    mkdir_p(paths.refsel_dir)

    # As an integrity check, we create a file with specific content (database version string) only after
    # everything else successfully completes.  Make sure this file does not exist until then.
    if os.path.isfile(paths.db_format):
        os.remove(paths.db_format)
    assert not os.path.isfile(paths.db_format)

    shutil.copyfile(refs_fasta, paths.refs_fasta)

    kmer_sizes = util.misc.make_seq(kmer_sizes or (27 if clark_variant=='light' else 31))
    dump_file(paths.n_segments_per_genome, n_segments_per_genome)
    dump_file(paths.kmer_sizes, ','.join(map(str, kmer_sizes)))

    mkdir_p(paths.clark_dir)

    # from CLARK manual:

    ##############
    # First, the user must define targets for the classification (the "targets_addresses.txt"
    # file). For metagenomics, for a classification at the species-level, then the user must
    # define targets accordingly.  So, at the species level, the targets_addresses.txt file is:
    #
    # $ cat targets_addresses.txt
    # /home/user1/bin/genome_1.fa S1
    # /home/user1/bin/genome_2.fa S1
    # /home/user1/bin/genome_3.fa S1
    # /home/user1/bin/genome_4.fa S2
    # /home/user1/bin/genome_5.fa S2
    # /home/user1/bin/genome_6.fa S2
    # /home/user1/bin/genome_7.fa S3
    # /home/user1/bin/genome_8.fa S3
    # /home/user1/bin/genome_9.fa S3
    # /home/user1/bin/genome_10.fa S3
    # $
    #
    # In this case, there are three distinct targets (S1, S2 and S3) that the program
    # will use for the classification.
    ##############

    # We're given the set of (possibly multi-segment) reference genomes as one fasta file.
    # We split it up into one fasta file per segment of a reference.
    # (CLARK needs single-sequence fastas as input).

    mkdir_p(paths.refs_segments_dir)

    with open(paths.clark_targets_addresses, 'wt') as out_targets_addresses:
        segments_all = tuple(Bio.SeqIO.parse(paths.refs_fasta, 'fasta'))
        assert (len(segments_all) % n_segments_per_genome) == 0, \
            "Number of reference sequences not a multiple of the number of reference segments"
        n_refs = int(len(segments_all) / n_segments_per_genome)

        _log.info('Processing %d references, each with %d segments', n_refs, n_segments_per_genome)

        for ref_num in range(n_refs):
            ref_idx_in_list = ref_num*n_segments_per_genome
            # we'll denote a ref by the id of its first segment
            ref_id = segments_all[ref_idx_in_list].id

            _log.info('Processing reference %s (%d/%d)', ref_id, ref_num, n_refs)

            # save each segment of its reference in a separate fasta
            # for CLARK, write out the correspondence between segment filenames and reference IDs
            # (in CLARK terminology each reference is a "target object")
            for segment_num in range(n_segments_per_genome):
                segment_rec = segments_all[ref_idx_in_list+segment_num]
                segment_id = segment_rec.id

                segment_fname = paths.ref_segment_fasta(ref_id, segment_id)
                Bio.SeqIO.write([segment_rec], segment_fname, 'fasta')

                clark_target_id = 'REF{:04d}_{}'.format(ref_num, ref_id)[:25]
                out_targets_addresses.write(segment_fname + ' ' + clark_target_id + '\n')

        # end: for ref_num in range(n_refs)

    dump_file(paths.clark_variant, clark_variant)

    clark = tools.clark.ClarkTool()

    for kmer_size in kmer_sizes:
        _log.info('Building CLARK database for kmer size %d', kmer_size)
        mkdir_p(paths.clark_db_dir(kmer_size))
        with util.file.tempfname(suffix='clark-db-bld') as dummy_result_fname:
            clark.execute(['-T', paths.clark_targets_addresses,
                           '-D', paths.clark_db_dir(kmer_size)+'/',
                           # CLARK does not have a separate "build-the-database" command; you give it a 
                           # FASTA file to classify, and then if the database isn't built yet it'll build it
                           # before classifying the file.  Here we just need the database, so we give CLARK a
                           # random single-segment file to classify, and discard the result.
                           '-O', segment_fname,
                           '-R', dummy_result_fname,
                           '-m', '0',
                           '-k', str(kmer_size)], variant=clark_variant)

    dump_file(paths.db_format, _db_format_string())

# end: def build_refsel_db(self, refs_fasta, n_segments_per_genome, kmer_sizes)

def _sanity_check_ref_genomes(refs_fasta, n_segments_per_genome, segment_len_margin=.1):
    '''Sanity checks: check for duplicate seqs or incomplete multi-segment genomes; check that within each
    segment number, all segments with that number have similar size (since viruses tend not to have
    long indels).

    Args:
        refs_fasta: fasta file of reference genomes, with segments of each multi-segment genome adjacent and
            with segments of each genome in the same segment order
        n_segments_per_genome: number of segments in the genome
        segment_len_margin: for each segment number, all reference segments with that segment number must have length
           within this fraction of the mean
    '''
    refsel_ids_seen = set()
    refsel_seqs_seen = set()
    segment_num_to_lens = collections.defaultdict(list)
    for rec_num, rec in enumerate(Bio.SeqIO.parse(refs_fasta, 'fasta')):
        assert rec.id not in refsel_ids_seen, 'Duplicate entry ID {} in {}'.format(rec.id, refs_fasta)
        refsel_ids_seen.add(rec.id)
        segment_seq = str(rec.seq)
        assert segment_seq not in refsel_seqs_seen, 'Duplicate seq under {} in {}'.format(rec.id, refs_fasta)
        segment_num_to_lens[rec_num % n_segments_per_genome].append(len(segment_seq))
    assert refsel_ids_seen, 'Empty reference list'
    assert (len(refsel_ids_seen) % n_segments_per_genome) == 0, \
        'Number of fasta records must be a multiple of n_segments_per_genome'
    for lens in segment_num_to_lens.values():
        mean_len = util.misc.mean(lens)
        min_len = mean_len*(1-segment_len_margin)
        max_len = mean_len*(1+segment_len_margin)
        assert all(min_len < L < max_len for L in lens), 'Inconsistent segment lengths'

# end: _sanity_check_ref_genomes(self, refs_fasta, n_segments_per_genome, segment_len_margin=.1)

def parser_build_refsel_db(parser=argparse.ArgumentParser()):
    '''Create parser for build_refsel_db command'''
    parser.add_argument('refs_fasta',
                        help='Reference genomes; segments of each genome must be consequtive fasta records, '
                        'always listed in the same segment order.')
    parser.add_argument('n_segments_per_genome', type=int, help='Number of segments in each reference genome.')
    parser.add_argument('refsel_dir', help='Directory where to put the refsel database.')
    parser.add_argument('--clark_variant', default='light', 
                        help='CLARK tool variant to use ("standard" or "light"); light variant, the default, '
                        'uses less memory but is less precise')
    util.cmd.attach_main(parser, build_refsel_db, split_args=True)
    return parser

__commands__.append(('build_refsel_db', parser_build_refsel_db))

# ====================
# ***  select_ref  ***
# ====================

def select_ref(refsel_dir, contigs, out_ref, min_contig_len=200):  # pylint: disable=R0914
    '''Choose a reference closest to the sample.

    Args:
        refsel_dir: directory where the refsel database has been previously stored by the build_refsel_db
            command.
        contigs: fasta file of contigs de novo assembled from the sample.  Repeated and overlapping contigs
            are ok; it is the set of kmers in the contigs that matters, not the frequency of each kmer.
        out_ref: the selected reference is copied here
        min_contig_len: contigs shorter than this are ignored
    '''

    from util.file import slurp_file

    paths = _RefSelPaths(refsel_dir)
    del refsel_dir  # make sure all paths are accessed via the `paths` object

    assert slurp_file(paths.db_format) == _db_format_string(), 'Database format mismatch'

    clark_variant = slurp_file(paths.clark_variant)

    # for now, just use the first kmer size
    kmer_size = int(slurp_file(paths.kmer_sizes).strip().split(',')[0])

    clark = tools.clark.ClarkTool()
    n_segments_per_genome = int(slurp_file(paths.n_segments_per_genome).strip())

    with util.file.tempfnames(('refsel_contigs_joined', 'refsel_clark_results')) as \
         (contigs_joined_fname, clark_results_fname):
        # compute a unified set of kmers from all contigs, by joining all contigs into one long string
        # (separated by runs of N kmer_size long to avoid creating chimeric kmers)
        contigs_joined = ('N'*kmer_size).join([str(r.seq) for r in Bio.SeqIO.parse(contigs, 'fasta') \
                                               if len(r.seq) > min_contig_len])
        Bio.SeqIO.write([Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(contigs_joined, Bio.Alphabet.generic_dna),
                                                 id='joined_contigs')], contigs_joined_fname, 'fasta')

        _log.info('running CLARK...')
        clark.execute(['-T', paths.clark_targets_addresses,
                       '-D', paths.clark_db_dir(kmer_size)+'/',
                       '-O', contigs_joined_fname,
                       '-m', '0', '--extended',
                       '-k', str(kmer_size),
                       '-R', clark_results_fname], variant=clark_variant)
        _log.info('CLARK done')
        closest_ref_contents = slurp_file(clark_results_fname+'.csv')
        closest_ref_lines = closest_ref_contents.split('\n')
        closest_ref_data_line = closest_ref_lines[1]
        closest_ref_name_list = closest_ref_data_line.split(',')
        closest_ref_name = closest_ref_name_list[-5]
        found = False
        for target_addr_lines in util.misc.grouper(slurp_file(paths.clark_targets_addresses).strip().split('\n'),
                                                   n_segments_per_genome):
            ref_fasta, ref_name = target_addr_lines[0].strip().split()
            if ref_name == closest_ref_name:
                util.file.concat([target_addr_line.strip().split()[0]
                                  for target_addr_line in target_addr_lines], out_ref)
                _log.info('Copying best ref %s to %s', ref_fasta, out_ref)
                found = True
                break
        assert found, 'Closest reference not found -- should not happen'

# end: def select_ref(refsel_dir, contigs, out_ref)

def parser_select_ref(parser=argparse.ArgumentParser(description='Select reference closest to sample')):
    '''Create parser for select_ref command'''
    parser.add_argument('refsel_dir', help='Directory where refsel info was stored by the build_refsel_db command.')
    parser.add_argument('contigs', help='Contigs assembled from sample.')
    parser.add_argument('out_ref', help='Output the selected reference to this file.')
    parser.add_argument('--min_contig_len', type=int, default=200, help='Ignore contigs shorter than this.')
    util.cmd.attach_main(parser, select_ref, split_args=True)
    return parser

__commands__.append(('select_ref', parser_select_ref))

def full_parser():
    '''Create parser for refsel commands'''
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
