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

    shutil.copyfile(refFasta, os.path.join(scafrefDir, 'reference.fasta'))


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
    kmc.dump(scafrefKmcDb, scafrefKmcDb+'.kmers.txt')
    for cx in (1,2,3):
        kmc.reduce(kmerDbIn=scafrefKmcDb, kmerDbOut=os.path.join(kmcDir, 'scafref.cx{}.kmc'.format(cx)),
                   kmerDbIn_opts='-ci1 -cx{}'.format(cx))

    clarkTargsFN=os.path.join(clarkDir,'targets_addresses.txt')

    segsRecs = [ [] for seg_idx in range(nsegs) ]

    refListFN=os.path.join(scafrefDir, 'reflist.txt')

    with open(clarkTargsFN, 'wt') as outTargAddr, open(refListFN, 'wt') as refListOut:
        recs = list(Bio.SeqIO.parse(scafrefFN, 'fasta'))
        assert (len(recs) % nsegs)==0, "Number of ref seqs not a multiple of the number of ref segments"
        nrefs = int(len(recs) / nsegs)

        log.info('Processing {} scaffolding references, with {} segments each'.format(nrefs,nsegs))

        for ref_num in range(nrefs):
            ref_idx_in_list = ref_num*nsegs
            # we'll denote a ref by the id of its first segment
            ref_id = recs[ref_idx_in_list].id

            log.info('Processing ref {} ({}/{})'.format(ref_id,ref_num,nrefs))
            
            refFNbase=util.file.string_to_file_name(ref_id)
            refFN=os.path.join(refsDir, refFNbase+'.fa')
            Bio.SeqIO.write(recs[ref_idx_in_list:ref_idx_in_list+nsegs], refFN, 'fasta')
            kmc.build_kmer_db(kmerDb=os.path.join(kmcDir, refFNbase+'.kmc'),
                              inFiles=refFN, kmerSize=kmcKmerSize,kmerOpts='-ci1')
            
            segFNs=[]
            for seg_idx in range(nsegs):
                seg_rec = recs[ref_idx_in_list+seg_idx]
                segsRecs[seg_idx].append(seg_rec)
                seg_id = seg_rec.id


                segFN=os.path.join(refsSplitDir,util.file.string_to_file_name(seg_id) + '.fa')
                segFNs.append(segFN)
                log.info('Procesing SEG {} ref_num={} ref_idx_in_list={} seg_id={} segFN={}'.format(seg_idx, ref_num, ref_idx_in_list, seg_id, segFN))
                Bio.SeqIO.write([seg_rec], segFN, 'fasta')
                targId='REF{:04d}_{}'.format(ref_num,refFNbase)
                targId=targId[:25]
                outTargAddr.write( segFN + ' ' + targId + '\n' )

            refListOut.write('\t'.join([refFN] + segFNs) + '\n')


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



def kmc_build_db(inFile, kmerDb, kmerSize, kmerOpts, histogram='', dump=''):
    '''Build a kmc kmer database'''

    kmc = tools.kmc.KmcTool()
    kmc.build_kmer_db(inFiles=[inFile], kmerDb=kmerDb, kmerSize=kmerSize, kmerOpts=kmerOpts)
    if histogram: kmc.histogram(kmerDb, histogram)
    if dump: kmc.dump(kmerDb, dump)
    
    
def parser_kmc_build_db(parser=argparse.ArgumentParser()):
    parser.add_argument('inFile', help='Bam or fasta file to go into the database.')
    parser.add_argument('kmerDb', help='kmer db name.')
    parser.add_argument('kmerSize', help='kmer size.')
    parser.add_argument('--kmerOpts', help='kmer opts.')
    parser.add_argument('--histogram', help='also output histogram.')
    parser.add_argument('--dump', help='also output dump.')
    util.cmd.attach_main(parser, kmc_build_db, split_args=True)
    return parser


__commands__.append(('kmc_build_db', parser_kmc_build_db))


def subtract_reads_from_refs(scafrefDir, readsKmcDb):
    '''For each ref, see which kmers in that ref are _not_ in the db'''
    
    


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
