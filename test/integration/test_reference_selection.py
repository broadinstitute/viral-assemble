'''Unit tests for automatic reference selection'''

__author__ = "ilya@broadinstitute.org"

import os.path
import tempfile

import Bio.SeqIO

import util.file
import reference_selection

def test_refsel_rbv():
    '''Test reference selection on the rabies virus'''

    from os.path import join
    
    work_dir = tempfile.mkdtemp()
    rbv_refsel_dir = join(work_dir, 'rbv_refsel_db')
    rbv_refs_fasta = join(work_dir, 'rbv_refs.fasta')
    test_dir=join(util.file.get_test_input_path(), 'TestReferenceSelection')

    util.file.gunzip(join(test_dir, 'all-rabies-genomes-t.fasta.gz'), rbv_refs_fasta)
    
    reference_selection.build_refsel_db(refs_fasta=rbv_refs_fasta, n_segments_per_genome=1, refsel_dir=rbv_refsel_dir,
                                        clark_variant='-l')
    rbv_selected_ref_fasta = join(work_dir, 'rbv_selected_ref.fasta')
    reference_selection.select_ref(refsel_dir=rbv_refsel_dir, contigs=join(test_dir, 'RBV51.assembly1-trinity.fasta'),
                                   out_ref=rbv_selected_ref_fasta)
    rbv_selected_ref = list(Bio.SeqIO.parse(rbv_selected_ref_fasta, 'fasta'))
    assert len(rbv_selected_ref) == 1
    assert rbv_selected_ref[0].id == 'JQ685920'

# end: def test_refsel_rbv()
