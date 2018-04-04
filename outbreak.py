#!/usr/bin/env python

"""Better tracing of who infected whom in an outbreak.

Commands in this module deal with aggregating data from multiple samples, to better understand what species are 
circulating in an outbreak, and who likely infected whom.
"""

__commands__ = []

# built-ins
import argparse
import logging
import glob
import itertools
import functools
import operator
import csv
import datetime
import os.path
import subprocess

# intra-module
import util.cmd
import util.file
import util.misc
import tools.blast
import tools.kmc
import read_utils

# third-party
import numpy
import Bio.SeqIO
import Bio.Blast
import Bio.Blast.Applications
import Bio.Blast.NCBIXML
import Bio.Blast.NCBIWWW

import Bio.Alphabet

_log = logging.getLogger(__name__)

# def gather_contigs_from_samples(sample_contig_fastas, out_fasta, min_contig_len=100):
#     """Given a list of per-sample fastas containing the contigs for each sample,
#     gather them into a single fasta, while changing the sequence ID line to indicate the
#     origin of each sequence.  Sample name is assumed to be the filename up to the first '.'.
#     Contigs shorter than `min_contig_len` are discarded."""


def gather_seq_names(in_fasta, out_seqnames):
    """Gather the sequence names"""
    seq_recs = tuple(Bio.SeqIO.parse(in_fasta, 'fasta'))
    print('got {} seq recs'.format(len(seq_recs)))
    sample_names = sorted([seq_rec.id[:seq_rec.id.index('|')] for seq_rec in seq_recs])
    print('got {} sample names'.format(len(sample_names)))
    util.file.dump_file(out_seqnames, '\n'.join(sample_names))


def gather_kmc_dbs():
    """Compute a kmc db for each sample"""

    samples = util.file.slurp_file('/idi/sabeti-scratch/ilya/sw/zvir/viral-ngs-etc/projects/mumps/tmp/seqnames.txt').strip().split()
    util.file.mkdir_p('tmp/kmc/sample_contig_kmers')
    kmcTool = tools.kmc.KmcTool()
    for sample in samples:
        kmcTool.build_db(['-k25', '-fm', '-cs1', 'tmp/sampnamed/{}.assembly1-spades.fasta'.format(sample),
                          'tmp/kmc/sample_contig_kmers/{}'.format(sample)])

def compute_kmc_union():
    """Compute the union, giving for each kmer in how many samples it is"""
    samples = util.file.slurp_file('/idi/sabeti-scratch/ilya/sw/zvir/viral-ngs-etc/projects/mumps/tmp/seqnames.txt').strip().split()
    kmcTool = tools.kmc.KmcTool()
    
    with open('tmp/kmc/unionfile.txt', 'w') as out:
        out.write('INPUT:\n')
        for sample in samples:
            sample_alpha = sample.replace('-', '_')
            out.write('{} = tmp/kmc/sample_contig_kmers/{} -ci1\n'.format(sample_alpha, sample))
        out.write('OUTPUT:\n')
        out.write('tmp/kmc/uniondb = {}\n'.format(' + '.join([sample.replace('-', '_') for sample in samples])))
        out.write('OUTPUT_PARAMS:\n')
        out.write('-ci1')
        out.write('\n')

def gather_assemblies():
    """Gather assemblies into one place"""
    samples = util.file.slurp_file('/idi/sabeti-scratch/ilya/sw/zvir/viral-ngs-etc/projects/mumps/tmp/seqnames.txt').strip().split()
    util.file.mkdir_p('tmp/sampnamed')
    util.file.mkdir_p('tmp/over300')
    print('samples=', samples)
    for sample in samples:

        print('looking at {}'.format(sample))

        if os.path.isfile('tmp/over300/{}.assembly1-spades.fasta'.format(sample)): continue
        if not os.path.isfile('tmp/assemblies/{}.raw.cleaned.assembly1-spades.fasta'.format(sample)):
            subprocess.check_call('dx download -a assemblies-de_novo/{}.raw.cleaned.assembly1-spades.fasta -o tmp/assemblies'.format(sample), shell=True)
        read_utils.prepend_to_fasta_headers('tmp/assemblies/{}.raw.cleaned.assembly1-spades.fasta'.format(sample),
                                            'lcl|'+sample+'|', 'tmp/sampnamed/{}.assembly1-spades.fasta'.format(sample))
        read_utils.filter_fasta_seqs('tmp/sampnamed/{}.assembly1-spades.fasta'.format(sample),
                                     'tmp/over300/{}.assembly1-spades.fasta'.format(sample),
                                     min_len=300)

    util.file.concat(['tmp/over300/{}.assembly1-spades.fasta'.format(sample) for sample in samples],
                     'tmp/all.over300.fasta')
    

def gather_infector_stats(seqs_fasta, seqs_dates, chainfiles, outfile, n_burnin, pair_reporting_threshold):
    """For each pair of samples (i,j), gather the probability that i infected j.

    Args:
        seqs_fasta: fasta file with the sequences (we only care about their names and count)
        seqs_dates: file with collection dates (just for joining to output)
        chainfiles: list of outbreaker chain files
        outfile: tsv file for output
        n_burnin: number of initial lines from each chain file to ignore
        pair_reporting_threshold: report only sample pairs that appear in this fraction of mcmc samples as infecting 
           one another
    """

    seqs = tuple(Bio.SeqIO.parse(seqs_fasta, 'fasta'))
    seqs_dates = util.file.slurp_file(seqs_dates).strip().split()
    assert len(seqs) == len(seqs_dates)

    n_samples = None
    n_mcmc_samples = 0
    for chainfile in chainfiles:
        with open(chainfile) as f:
            headers = f.readline().strip().split()
            alpha_cols = [i for i, h in enumerate(headers) if h.startswith('alpha_')]
            if n_samples is None:
                n_samples = len(alpha_cols)
                assert n_samples == len(seqs)
                infector_counts = numpy.zeros((n_samples, n_samples), dtype=int)
            else:
                assert len(alpha_cols) == n_samples

        assert headers[alpha_cols[0]] == 'alpha_1'
        assert headers[alpha_cols[-1]] == 'alpha_{}'.format(n_samples)
        t = numpy.loadtxt(chainfile, dtype=int, skiprows=1+n_burnin, usecols = alpha_cols)
        #numpy.savetxt('/tmp/f.txt', t, fmt='%d', delimiter='\t')
        n_mcmc_samples += t.shape[0]
        print(type(t), t.dtype, t.shape)
        for sample in range(n_samples):
            for infector, count in enumerate(numpy.bincount(t[:,sample])):
                if infector > 0:
                    infector_counts[sample, infector-1] += count

    # end: for chainfile in chainfiles

    #numpy.savetxt('/tmp/f.txt', infector_counts, fmt='%d',
    #              header='\t'.join('s{}'.format(s+1) for s in range(n_samples)), delimiter='\t')

    n_mcmc_threshold = int(n_mcmc_samples * pair_reporting_threshold)
    print('n_mcmc_samples=', n_mcmc_samples, 'n_mcmc_threshold=', n_mcmc_threshold)
    
    assert (numpy.diagonal(infector_counts) == 0).all()

    pairs = []

    for i in range(n_samples):
        for j in range(i, n_samples):
            cnt_ij = infector_counts[i,j] + infector_counts[j,i]
            if cnt_ij > n_mcmc_threshold:
                pairs.append((float(cnt_ij) / n_mcmc_samples, (i,j)))


    def days_between(d1, d2):
        strptime = datetime.datetime.strptime
        d1 = strptime(d1, "%Y-%m-%d")
        d2 = strptime(d2, "%Y-%m-%d")
        return abs((d2 - d1).days)

    pair2info = {}

    with open(outfile, 'w') as out:
        w = csv.DictWriter(out, fieldnames=['prob', 'id1', 'id2', 'days_apart'], delimiter='\t')
        w.writeheader()
        for prob, (id1, id2) in sorted(pairs, reverse=True):
            days_apart = days_between(seqs_dates[id1], seqs_dates[id2])
            w.writerow({'prob': '{:.6f}'.format(prob), 'id1': seqs[id1].name, 'id2': seqs[id2].name,
                        'days_apart': days_apart})

            def get_sample(id): return id[:id.index('|')]
            id1_sample = get_sample(seqs[id1].name)
            id2_sample = get_sample(seqs[id2].name)
            pair2info[(id1_sample, id2_sample)] = {'prob': prob, 'days_apart': days_apart}
            pair2info[(id2_sample, id1_sample)] = {'prob': prob, 'days_apart': days_apart}

    return pair2info

    #print '\n'.join('{:.4f} {}'.format(cnt, pair) for cnt, pair in sorted(pairs))

def find_common_contigs(sample_contigs_fastas, min_contig_len=1000):
    """Given de novo assemblies of contigs from a set of samples, find contig parts that appear in 
    multiple samples.

    Implementation:

    Make a blast database from all contigs, then blast each contig against it.
    """

    # for each fasta, prepend a sample ID to the sequences in it.

    with util.file.tmp_dir(suffix='common_contigs_dir') as tdir:
        tdir = '/broad/hptmp/ilya'
        seqs = []

        sample_names = []
        sample_seqids = []
        sampnamed_fastas = []

        for sample_contigs_fasta in sample_contigs_fastas:
            contigs = tuple(Bio.SeqIO.parse(sample_contigs_fasta, 'fasta'))

            contigs = [s for s in contigs if len(s.seq) > min_contig_len]
        
            contig_seqids = []

            bname = os.path.basename(sample_contigs_fasta)
            sample_name = bname[:bname.index('.')]
            sample_names.append(sample_name)

            for seq_rec in contigs:
                seq_rec.id = sample_name + '.' + seq_rec.id
                seq_rec.name = sample_name + '.' + seq_rec.name
                seq_rec.description = sample_name + '.' + seq_rec.description
                contig_seqids.append(seq_rec.id)
            sample_seqids.append(contig_seqids)
                
            sampnamed_fastas.append(os.path.splitext(sample_contigs_fasta)[0] + '.sampnamed.fasta')
            print(sampnamed_fastas)
            Bio.SeqIO.write(contigs, sampnamed_fastas[-1], 'fasta')
            seqs.extend(contigs)
        seqs = [s for s in seqs if len(s.seq) > min_contig_len]
        seqs_fasta = os.path.join(tdir, 'seqs.fasta')
        Bio.SeqIO.write(seqs, seqs_fasta, 'fasta')
        
        db = os.path.join(tdir, 'seqs_blastdb')
        tools.blast.MakeblastdbTool().build_database([seqs_fasta], db, parse_seqids=True)

        blast = tools.blast.BlastnTool()
        for i, sample_contigs_fasta in enumerate(sample_contigs_fastas[:-1]):
            seqid_list_fname = os.path.join(tdir, 'seqidlist_{}.txt'.format(i))
            util.file.dump_file(seqid_list_fname, '\n'.join(functools.reduce(operator.concat,
                                                                             sample_seqids[i+1:], [])))
            blastn_cline = Bio.Blast.Applications.NcbiblastnCommandline(query=sampnamed_fastas[i], 
                                                                        db=db, evalue=0.000001,
                                                                        seqidlist=seqid_list_fname,
                                                                        ungapped=True,
                                                                        qcov_hsp_perc=20,
                                                                        outfmt=5,
                                                                        out=os.path.join(tdir, '{:03d}.{}.blast.xml'.format(i, sample_names[i])))
            print(blastn_cline)
            blastn_cline()



def get_top_hits(pair2info):
    """Find top blast hits"""

    tdir = '/broad/hptmp/ilya'
    for xmlfname in sorted(set(glob.glob('/broad/hptmp/ilya/???.MuV-???.blast.xml'))):
        with open(xmlfname) as xml_f:
            blast_recs = list(Bio.Blast.NCBIXML.parse(xml_f))
            print(xmlfname, len(blast_recs))
            for blast_rec in blast_recs:
                id1 = blast_rec.query[:20]
                for alignment in blast_rec.alignments:
                    id2 = alignment.hit_id[:20]
                    hsps = [hsp.align_length for hsp in alignment.hsps if hsp.bits > 4000 and hsp.align_length > 500]
                    if hsps:
                        print(id1, id2, hsps, pair2info.get((id1[:7], id2[:7]), None))


                      
            
    
#    with util.file.tmp_dir(suffix='contigs_blast_db') as blastdb_dir:
#        tools.blast.MakeblastdbTool().build_database(


if __name__ == '__main__':
    if False:
        chains_dir = '/idi/sabeti-scratch/swohl/mumps/outbreaker/v3-final-SH'
        chains_files = glob.glob(os.path.join(chains_dir, 'run?-II-outbreak.SH.rep?.1e6.chains.txt'))
        seqs_fasta = os.path.join(chains_dir, 'II-outbreak-v3.SH.aligned.pruned.fasta')
        seqs_dates = os.path.join(chains_dir, 'II-outbreak-v3.dates.txt')
        print(chains_files)
        pair2info = gather_infector_stats(seqs_fasta, seqs_dates, chains_files, 'infection_pairs.tsv', 
                                          n_burnin=51, pair_reporting_threshold=.01)


    if False:
        samples = []
        with open('infection_pairs.tsv') as top_pairs_f:
            reader = csv.DictReader(top_pairs_f, delimiter='\t')
            for row in reader:
                if int(row['days_apart']) < 35:
                    samples.append(row['id1'][:7])
                    samples.append(row['id2'][:7])
        samples = sorted(list(set(samples)))

        contigs_dir = '/idi/sabeti-scratch/ilya/sw/zvir/viral-ngs-etc/projects/mumps/tmp/assemblies'
        #contigs_fastas = sorted(glob.glob(os.path.join(contigs_dir, 'MuV-*.raw.cleaned.assembly1-spades.fasta')))
        contigs_fastas = sorted([os.path.join(contigs_dir, 
                                              '{}.raw.cleaned.assembly1-spades.fasta'.format(s)) for s in samples])
        contigs_fastas = list(filter(os.path.isfile, contigs_fastas))

        print('\n'.join(contigs_fastas))
        find_common_contigs(contigs_fastas)

    #get_top_hits(pair2info)
    #gather_seq_names(in_fasta='/idi/sabeti-scratch/swohl/mumps/outbreaker/v3-final/II-outbreak-v3.aligned.pruned.fasta',
    #                 out_seqnames='tmp/seqnames.txt')
    #gather_assemblies()
    #gather_kmc_dbs()
    compute_kmc_union()





