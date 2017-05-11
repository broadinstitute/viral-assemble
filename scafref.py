#!/usr/bin/env python
''' Management of scaffolding references (scafrefs).
'''

__author__ = "ilya@broadinstitute.org"
__commands__ = []

# built-ins
import argparse
import logging
import random
import os
import os.path
import re
import shutil
import subprocess

try:
    from itertools import zip_longest    # pylint: disable=E0611
except ImportError:
    from itertools import izip_longest as zip_longest    # pylint: disable=E0611

# intra-module
import util.cmd
import util.file
import util.genbank
import util.misc
import util.vcf
import read_utils
import tools
import tools.picard
import tools.samtools
import tools.gatk
import tools.novoalign
import tools.trimmomatic
import tools.trinity
import tools.mafft
import tools.mummer
import tools.muscle
import tools.spades
import tools.kmc
import tools.clark

# third-party
import Bio.AlignIO
import Bio.SeqIO
import Bio.Data.IUPACData

log = logging.getLogger(__name__)

def gather_refs(refFasta, lastalFasta, refsFasta):
    '''Gather possible reference genomes on which to scaffold'''
    util.file.concat([refFasta,lastalFasta], refsFasta)


def parser_gather_refs(parser=argparse.ArgumentParser()):
    parser.add_argument('refFasta', help='Main ref fasta.')
    parser.add_argument('lastalFasta', help='Lastal seqs.')
    parser.add_argument('refsFasta', help='Output of all scaffolding refs.')
    util.cmd.attach_main(parser, gather_refs, split_args=True)
    return parser


__commands__.append(('gather_refs', parser_gather_refs))



def build_scafref_db(scafrefDir, refFasta, kmcKmerSize=35, clarkKmerSize=35,
                     JVMmemory=None,
                     threads=1,
                     novoalign_license_path=None
                 ):
    '''Precompute info that will help us identify the reference to align to.

    We start with the file {scafrefDir}/scafref.fasta, which contains the possible references to which
    we can scaffold.  The reference consists of nseg segments.  scafref.fasta looks like

    ref1_seg1
    ref1_seg2
    ...
    ref1_segN
    ref2_seg1
    ref2_seg2
    ...
    ref2_segN
    '''

    scafrefFN=os.path.join(scafrefDir,'scafref.fasta')

    doneFile=scafrefFN+'.indexing.done'
    if os.path.isfile(doneFile): os.unlink(doneFile)

    nsegs=util.file.fasta_length(refFasta)

    read_utils.index_fasta_all(scafrefFN, JVMmemory=JVMmemory, threads=threads, 
                               novoalign_license_path=novoalign_license_path)

    clarkDir=os.path.join(scafrefDir,'CLARK', 'k'+str(clarkKmerSize))
    clarkDbDir=os.path.join(clarkDir,'clarkDb')
    util.file.mkdir_p(clarkDbDir)

    refsDir=os.path.join(scafrefDir,'refs')
    util.file.mkdir_p(refsDir)
    refsSplitDir=os.path.join(clarkDir,'refs')
    util.file.mkdir_p(refsSplitDir)
    kmcDir=os.path.join(scafrefDir,'kmc','k'+str(kmcKmerSize))
    util.file.mkdir_p(kmcDir)

    kmc = tools.kmc.KmcTool()

    scafrefKmcDb=os.path.join(kmcDir, 'scafref.kmc')

    kmc.build_kmer_db(kmerDb=scafrefKmcDb,
                      inFiles=scafrefFN, kmerSize=kmcKmerSize,kmerOpts='-ci1')
    kmc.histogram(scafrefKmcDb, scafrefKmcDb+'.histogram.txt')
    kmc.histogram(scafrefKmcDb, scafrefKmcDb+'.kmers.txt')
    for cx in (1,2,3):
        kmc.execute(tool_sfx='_tools',
                    args=['reduce',os.path.join(kmcDir,'scafref'.format(cx)),
                          '-ci1','-cx{}'.format(cx), 
                          os.path.join(kmcDir,'scafref.cx{}'.format(cx))])

    clarkTargsFN=os.path.join(clarkDir,'targets_addresses.txt')

    segsRecs = [ [] for seg_idx in range(nsegs) ]


    with open(clarkTargsFN, 'wt') as outTargAddr:
        recs = list(Bio.SeqIO.parse(scafrefFN, 'fasta'))
        assert (len(recs) % nsegs)==0, "Number of ref seqs not a multiple of the number of ref segments"
        nrefs = int(len(recs) / nsegs)

        log.info('Processing {} scaffolding references, with {} segments each'.format(nrefs,nsegs))

        for ref_num in range(nrefs):
            ref_idx_in_list = ref_num*nsegs
            # we'll denote a ref by the id of its first segment
            ref_id = util.genbank.parse_accession_str(recs[ref_idx_in_list].id)

            log.info('Processing ref {} ({}/{})'.format(ref_id,ref_num,nrefs))

            refFN=os.path.join(refsDir, util.file.string_to_file_name(ref_id)+'.fa')
            Bio.SeqIO.write(recs[ref_idx_in_list:ref_idx_in_list+nsegs], refFN, 'fasta')
            kmc.build_kmer_db(kmerDb=os.path.join(kmcDir, util.file.string_to_file_name(ref_id)+'.kmc'),
                              inFiles=refFN, kmerSize=kmcKmerSize,kmerOpts='-ci1')

            for seg_idx in range(nsegs):
                seg_rec = recs[ref_idx_in_list+seg_idx]
                segsRecs[seg_idx].append(seg_rec)
                seg_id = util.genbank.parse_accession_str(seg_rec.id)
                segFN=os.path.join(refsSplitDir,util.file.string_to_file_name(seg_id) + '.fa')
                Bio.SeqIO.write([seg_rec], segFN, 'fasta')
                targId='REF{:04d}_{}'.format(ref_num,util.file.string_to_file_name(seg_id))
                targId=targId[:25]
                outTargAddr.write( segFN + ' ' + targId + '\n' )

    for seg_idx in range(nsegs):
        segFasta=os.path.join(scafrefDir,'segs{}.fasta'.format(seg_idx))
        Bio.SeqIO.write(segsRecs[seg_idx], segFasta, 'fasta')
        read_utils.index_fasta_all(segFasta, JVMmemory=JVMmemory, threads=threads, 
                                   novoalign_license_path=novoalign_license_path)

    clark = tools.clark.ClarkTool()

    clark.execute(['-T', clarkTargsFN,
                   '-D', clarkDbDir+'/',
                   '-O', os.path.join(scafrefDir, 'scafref.fasta'),
                   '-m', '0', '--extended',
                   '-R', os.path.join(clarkDir, 'scafref.test')])

    shutil.copyfile(scafrefFN, doneFile)
    
    # kmc = tools.kmc.KmcTool()
    
    # kmc.builds_kmer_db( inFiles = [refsFasta], 
    #                     kmerDb = kmerDb, kmerSize=kmerSize, kmerOpts=kmerOpts, threads=threads )

    
def parser_build_scafref_db(parser=argparse.ArgumentParser()):
    parser.add_argument('scafrefDir', help='Directory where scafref info is stored.')
    parser.add_argument('refFasta', help='Main fasta reference.')
    util.cmd.attach_main(parser, build_scafref_db, split_args=True)
    return parser


__commands__.append(('build_scafref_db', parser_build_scafref_db))


def choose_scafref(scafrefDir, ):
    '''Choose a reference to which to scaffold ccontigs'''
    
    clark = tools.clark.ClarkTool()

    make_clark_db(refsFasta=os.path.join(scafrefDir, 'scafref.fasta'), 
                  clarkDbDir=scafrefDir)

    clark.execute(['-T', os.path.join(scafrefDir, 'targets_addresses.txt'),
                   '-D', os.path.join(scafrefDir, 'clarkInternalDB/'),
                   '-O', os.path.join(scafrefDir, 'scafref.fasta'),
                   '-m', '0', '--extended',
                   '-R', os.path.join(scafrefDir, 'scafref.test')])
    
    
def parser_choose_scafref(parser=argparse.ArgumentParser()):
    parser.add_argument('scafrefDir', help='Directory where scafref info is stored.')
    util.cmd.attach_main(parser, choose_scafref, split_args=True)
    return parser


__commands__.append(('choose_scafref', parser_choose_scafref))


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
