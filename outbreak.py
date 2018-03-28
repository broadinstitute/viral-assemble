#!/usr/bin/env python

"""Better tracing of who infected whom in an outbreak"""

__commands__ = []

# built-ins
import argparse
import logging
import glob
import itertools
import csv
import datetime
import os.path

# intra-module
import util.cmd
import util.file
import util.misc
import tools.blast

# third-party
import Bio.SeqIO
import numpy

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

    print '\n'.join('{:.4f} {}'.format(cnt, pair) for cnt, pair in sorted(pairs))

def find_common_contigs(sample_contig_fastas):
    """Given de novo assemblies of contigs from a set of samples, find contig parts that appear in 
    multiple samples.

    Implementation:

    Make a blast database from all contigs, then blast each contig against it.
    """
    
#    with util.file.tmp_dir(suffix='contigs_blast_db') as blastdb_dir:
#        tools.blast.MakeblastdbTool().build_database(


if __name__ == '__main__':
    chains_dir = '/idi/sabeti-scratch/swohl/mumps/outbreaker/v3-final-SH'
    chains_files = glob.glob(os.path.join(chains_dir, 'run?-II-outbreak.SH.rep?.1e6.chains.txt'))
    seqs_fasta = os.path.join(chains_dir, 'II-outbreak-v3.SH.aligned.pruned.fasta')
    seqs_dates = os.path.join(chains_dir, 'II-outbreak-v3.dates.txt')
    print(chains_files)
    gather_infector_stats(seqs_fasta, seqs_dates, chains_files, 'infection_pairs.tsv', 
                          n_burnin=51, pair_reporting_threshold=.2)
