'''Tests of automatic reference selection'''

__author__ = "ilya@broadinstitute.org"

import logging

import Bio.SeqIO

import util.file
import reference_selection

_log = logging.getLogger(__name__)   # pylint: disable=C0103

def test_refsel(tmpdir):
    '''Test reference selection'''

    tests = {
        'rbv': {
            'n_segments_per_genome': 1,
            'samples':
            [{'contigs': 'RBV51.assembly1-trinity.fasta', 'best_ref': ['JQ685920']}]
        },
        'lasv': {
            'n_segments_per_genome': 2,
            'samples':
            [{'contigs': 'A4.assembly1-trinity.fasta',
              'best_ref': ['LASV-G502-S-Sierra_Leone-2009H', 'LASV-G502-L-Sierra_Leone-2009H']}]
        }
    }

    from os.path import join

    test_dir = join(util.file.get_test_input_path(), 'TestReferenceSelection')
    for vir_name, vir_info in sorted(tests.items()):
        _log.info('testing vir_name=%s', vir_name)

        work_dir = str(tmpdir.mkdir(vir_name))

        refsel_dir = join(work_dir, 'refsel_db')
        refs_fasta = join(work_dir, 'refs.fasta')

        util.file.gunzip(join(test_dir, vir_name, 'refs.fasta.gz'), refs_fasta)
        n_segments_per_genome = vir_info['n_segments_per_genome']
        reference_selection.build_refsel_db(refs_fasta=refs_fasta, n_segments_per_genome=n_segments_per_genome,
                                            refsel_dir=refsel_dir)
        selected_ref_fasta = join(work_dir, 'selected_ref.fasta')
        for sample in vir_info['samples']:
            reference_selection.select_ref(refsel_dir=refsel_dir,
                                           contigs=join(test_dir, vir_name, sample['contigs']),
                                           out_ref=selected_ref_fasta)
            selected_ref = list(Bio.SeqIO.parse(selected_ref_fasta, 'fasta'))
            assert len(selected_ref) == n_segments_per_genome
            assert [r.id for r in selected_ref] == sample['best_ref']

# end: def test_refsel()
