#!/usr/bin/env python

"""Better tracing of who infected whom in an outbreak"""

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

# intra-module
import util.cmd
import util.file
import util.misc
import tools.blast

# third-party
import numpy
import Bio.SeqIO
import Bio.Blast
import Bio.Blast.Applications

_log = logging.getLogger(__name__)

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

    with open(outfile, 'w') as out:
        w = csv.DictWriter(out, fieldnames=['prob', 'id1', 'id2', 'days_apart'], delimiter='\t')
        w.writeheader()
        for prob, (id1, id2) in sorted(pairs, reverse=True):
            w.writerow({'prob': '{:.4f}'.format(prob), 'id1': seqs[id1].name, 'id2': seqs[id2].name,
                        'days_apart': days_between(seqs_dates[id1], seqs_dates[id2])})

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
            
    
#    with util.file.tmp_dir(suffix='contigs_blast_db') as blastdb_dir:
#        tools.blast.MakeblastdbTool().build_database(


if __name__ == '__main__':
    if False:
        chains_dir = '/idi/sabeti-scratch/swohl/mumps/outbreaker/v3-final-SH'
        chains_files = glob.glob(os.path.join(chains_dir, 'run?-II-outbreak.SH.rep?.1e6.chains.txt'))
        seqs_fasta = os.path.join(chains_dir, 'II-outbreak-v3.SH.aligned.pruned.fasta')
        seqs_dates = os.path.join(chains_dir, 'II-outbreak-v3.dates.txt')
        print(chains_files)
        gather_infector_stats(seqs_fasta, seqs_dates, chains_files, 'infection_pairs.tsv', 
                              n_burnin=51, pair_reporting_threshold=.2)


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
