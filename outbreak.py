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
import collections

# intra-module
import util.cmd
import util.file
import util.misc
import tools.blast
import tools.kmc
import tools.cdhit
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
        kmcTool.build_kmer_db(['-k45', '-fm', '-cs1', 
                               'tmp/sampnamed/{}.assembly1-spades.fasta'.format(sample),
                               'tmp/kmc/sample_contig_kmers/k45.{}'.format(sample)])

def compute_kmc_union():
    """Compute the union, giving for each kmer in how many samples it is"""
    samples = util.file.slurp_file('/idi/sabeti-scratch/ilya/sw/zvir/viral-ngs-etc/projects/mumps/tmp/seqnames.txt').strip().split()
    kmcTool = tools.kmc.KmcTool()
    
    with open('tmp/kmc/unionfile.k45.txt', 'w') as out:
        out.write('INPUT:\n')
        for sample in samples:
            sample_alpha = sample.replace('-', '_')
            out.write('{} = tmp/kmc/sample_contig_kmers/k45.{} -ci1\n'.format(sample_alpha, sample))
        out.write('OUTPUT:\n')
        out.write('tmp/kmc/uniondb.k45 = {}\n'.format(' + '.join([sample.replace('-', '_') for sample in samples])))
        out.write('OUTPUT_PARAMS:\n')
        out.write('-ci1')
        out.write('\n')

def gather_kmc_dbs():
    """Compute a kmc db for each sample"""

    samples = util.file.slurp_file('/idi/sabeti-scratch/ilya/sw/zvir/viral-ngs-etc/projects/mumps/tmp/seqnames.txt').strip().split()
    util.file.mkdir_p('tmp/kmc/sample_contig_kmers')
    kmcTool = tools.kmc.KmcTool()
    for sample in samples:
        kmcTool.build_kmer_db(['-k45', '-fm', '-cs1', 
                               'tmp/sampnamed/{}.assembly1-spades.fasta'.format(sample),
                               'tmp/kmc/sample_contig_kmers/k45.{}'.format(sample)])


def find_linking_kmers_in_kmc_dbs():
    """For each sample, find kmers that link it with one or two others."""

    samples = util.file.slurp_file('/idi/sabeti-scratch/ilya/sw/zvir/viral-ngs-etc/projects/mumps/tmp/seqnames.txt').strip().split()
    util.file.mkdir_p('tmp/kmc/sample_contig_kmers')
    kmcTool = tools.kmc.KmcTool()
    for sample in samples:
        kmcTool.execute(['-v', 'simple',
                         'tmp/kmc/sample_contig_kmers/k45.{}'.format(sample), '-ci1',
                         'tmp/kmc/uniondb.k45', '-ci2', '-cx2', 'intersect',
                         'tmp/kmc/sample_contig_kmers/k45.in2samp.{}'.format(sample), '-ci1' ])
        kmcTool.execute(['-v', 'transform',
                         'tmp/kmc/sample_contig_kmers/k45.in2samp.{}'.format(sample), '-ci1',
                         'dump', 'tmp/kmc/sample_contig_kmers/k45.in2samp.{}.txt'.format(sample)])

def gather_linkingkmer_stats():
    """Gather stats about linking kmers"""
    samples = util.file.slurp_file('/idi/sabeti-scratch/ilya/sw/zvir/viral-ngs-etc/projects/mumps/tmp/seqnames.txt').strip().split()
    util.file.mkdir_p('tmp/kmc/sample_contig_kmers')
    kmcTool = tools.kmc.KmcTool()
    kmer2samps = collections.defaultdict(list)
    for sample in samples:
        lines = util.file.slurp_file('tmp/kmc/sample_contig_kmers/k45.in2samp.{}.txt'.format(sample)).strip().split('\n')
        kmers = [line.strip().split()[0] for line in lines]
        for kmer in kmers:
            assert len(kmer) == 45, 'kmer is {}'.format(kmer)
            kmer2samps[kmer].append(sample)
        print(sample, len(kmer2samps))
    print('got info for {} kmers'.format(len(kmer2samps)))
    pair2count = collections.Counter()
    for kmer, samps in kmer2samps.items():
        assert len(samps) == 2, 'kmer={} samps={}'.format(kmer, samps)
        pair2count[tuple(sorted(samps))] += 1
    spairs = sorted(pair2count.items(), key=operator.itemgetter(1), reverse=True)
    pair_counts_hist = collections.Counter()
    print('histogram=', numpy.histogram(list(map(operator.itemgetter(1), spairs))))
    with open('spairs.tsv', 'w') as out:
        out.write('id1\tid2\tnkmers\n')
        for pair, cnt in spairs:
            out.write('{}\t{}\t{}\n'.format(pair[0], pair[1], cnt))


        
################################################################################################
        
def gather_assemblies():
    """Gather assemblies into one place"""
    samples = util.file.slurp_file('/idi/sabeti-scratch/ilya/sw/zvir/viral-ngs-etc/projects/mumps/tmp/seqnames.txt').strip().split()
    util.file.mkdir_p('tmp/sampnamed')
    util.file.mkdir_p('tmp/over300')
    print('samples=', samples)
    for sample in samples:

        print('looking at {}'.format(sample))

#        if os.path.isfile('tmp/over300/{}.assembly1-spades.fasta'.format(sample)): continue
        if not os.path.isfile('tmp/assemblies/{}.raw.cleaned.assembly1-spades.fasta'.format(sample)):
            subprocess.check_call('dx download -a assemblies-de_novo/{}.raw.cleaned.assembly1-spades.fasta -o tmp/assemblies'.format(sample), shell=True)
        read_utils.prepend_to_fasta_headers('tmp/assemblies/{}.raw.cleaned.assembly1-spades.fasta'.format(sample),
                                            sample+'.', 'tmp/sampnamed/{}.assembly1-spades.fasta'.format(sample))
        read_utils.filter_fasta_seqs('tmp/sampnamed/{}.assembly1-spades.fasta'.format(sample),
                                     'tmp/over300/{}.assembly1-spades.fasta'.format(sample),
                                     min_len=300)

    util.file.concat(['tmp/over300/{}.assembly1-spades.fasta'.format(sample) for sample in samples],
                     'tmp/all.over300.fasta')
    

def gather_infector_stats(seqs_fasta, seqs_dates, chainfiles, outfile, n_burnin, pair_reporting_threshold):
    """For each pair of samples (i,j), gather the probability that i infected j, from the output of R package 'outbreaker'.

    Args:
        seqs_fasta: fasta file with the sequences (we only care about their names and the total number of seqs)
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


def make_samples_blastdbs():
    """Make a blast db for each sample"""
    samples = util.file.slurp_file('tmp/seqnames.txt').strip().split()
    util.file.mkdir_p('tmp/sampblast')
    for sample in samples:
        print(sample)
        tools.blast.MakeblastdbTool().build_database(['tmp/over300/{}.assembly1-spades.fasta'.format(sample)], 'tmp/sampblast/{}.blastdb'.format(sample), parse_seqids=True)

def match_contigs_from_two_samples(sample1, sample2):
    """Find close matches between contigs from two different samples.
    """
    samples = util.file.slurp_file('tmp/seqnames.txt').strip().split()
    sample1_fasta = 'tmp/sampnamed/{}.assembly1-spades.fasta'.format(sample1)
    sample2_fasta = 'tmp/sampnamed/{}.assembly1-spades.fasta'.format(sample2)

    util.file.mkdir_p('tmp/sampblast')

    cdhit = tools.cdhit.CdHit()
    tools.blast.MakeblastdbTool().build_database([contigs2+'_blastdb'], db, parse_seqids=True)
    blastn_cline = Bio.Blast.Applications.NcbiblastnCommandline(query=sampnamed_fastas[i], 
                                                                db=db, evalue=0.000001,
                                                                seqidlist=seqid_list_fname,
                                                                ungapped=True,
                                                                qcov_hsp_perc=20,
                                                                outfmt=5,
                                                                out=os.path.join(tdir, '{:03d}.{}.blast.xml'.format(i, sample_names[i])))    
    

def cmp_links():
    """Compare samples with links"""
    with open('/idi/sabeti-scratch/ilya/mlinks.csv') as mlinks:
        reader = csv.DictReader(mlinks)
        for row in reader:
            pass

def remove_mumps_contigs():
    """Remove mumps contigs (since we're asking whether other things get transferred along with mumps)."""


    with open('tmp/mumpscontigs.tsv') as mumpscontigs, open('tmp/mumpscontigs.filt.tsv', 'w') as mumpscontigsfilt:
        fieldnames='qseqid sseqid evalue bitscore length pident nident'.split()
        reader = csv.DictReader(mumpscontigs, delimiter='\t', fieldnames=fieldnames)
        writer = csv.DictWriter(mumpscontigsfilt, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        mumps_contigs = set()
        for row in reader:
            if row['qseqid'] in mumps_contigs: continue
            qseqid = row['qseqid'].split('_')
            contig_len = float(int(qseqid[3]))
            align_len = float(int(row['length']))
            qfrac = align_len / contig_len
            #assert qfrac <= 1
            if int(row['nident']) > 1000 and qfrac > .2:
                writer.writerow(row)
                mumps_contigs.add(row['qseqid'])

    samples = util.file.slurp_file('tmp/seqnames.txt').strip().split()
    util.file.mkdir_p('tmp/nomumps')

    for sample in samples:
        seq_recs = tuple(Bio.SeqIO.parse('tmp/over300/{}.assembly1-spades.fasta'.format(sample), 'fasta'))
        seq_recs_nomumps = [r for r in seq_recs if r.id not in mumps_contigs]
        if (len(seq_recs)-len(seq_recs_nomumps)) == 0:
            print('sample', sample, 'drop', len(seq_recs)-len(seq_recs_nomumps))
    
#    with util.file.tmp_dir(suffix='contigs_blast_db') as blastdb_dir:
#        tools.blast.MakeblastdbTool().build_database(

def plot_contact_tracing_info():
    """Plot the contact tracing info of connections between samples"""
    
    samples = util.file.slurp_file('/idi/sabeti-scratch/ilya/sw/zvir/viral-ngs-etc/projects/mumps/tmp/seqnames.txt').strip().split()

    with open('tmp/mlinks.csv') as mlinks, open('tmp/mlinks.dot', 'w') as dot:
        reader = csv.DictReader(mlinks, delimiter=',')
        dot.write('graph G {\n')
        def to_node(s): return s.replace('-', '_')
        for row in reader:
            if not (row['sample1'] in samples and row['sample2'] in samples): continue
            dot.write('{} -- {} [style={}];\n'.format(to_node(row['sample1']), to_node(row['sample2']),
                                                      'solid' if row['type']=='contact' else 'dashed'))
        dot.write('}\n')

    subprocess.check_call('dot -Tsvg -otmp/mlinks.svg tmp/mlinks.dot', shell=True)

def compute_num_diffs():
    """Compute diff counts"""
    seq_recs = tuple(Bio.SeqIO.parse('/idi/sabeti-scratch/swohl/mumps/outbreaker/v3-final/II-outbreak-v3.aligned.pruned.fasta', 'fasta'))
    muv_010 = [r for r in seq_recs if r.id.startswith('MuV-103')]
    muv_031 = [r for r in seq_recs if r.id.startswith('MuV-113')]
    assert len(muv_010) == len(muv_031) == 1
    assert len(muv_010[0].seq) == len(muv_031[0].seq)
    for i in range(len(muv_010[0].seq)):
        if muv_010[0].seq[i] != muv_031[0].seq[i]:
            print(i)


if __name__ == '__main__':
    if False:
        chains_dir = '/idi/sabeti-scratch/swohl/mumps/outbreaker/v3-final-SH'
        chains_files = glob.glob(os.path.join(chains_dir, 'run?-II-outbreak.SH.rep?.1e6.chains.txt'))
        seqs_fasta = os.path.join(chains_dir, 'II-outbreak-v3.SH.aligned.pruned.fasta')
        seqs_dates = os.path.join(chains_dir, 'II-outbreak-v3.dates.txt')
        print(chains_files)
        pair2info = gather_infector_stats(seqs_fasta, seqs_dates, chains_files, 'infection_pairs.SH.tsv', 
                                          n_burnin=101, pair_reporting_threshold=.7)

    if  False:
        chains_dir = '/idi/sabeti-scratch/swohl/mumps/outbreaker/v3-final'
        chains_files = glob.glob(os.path.join(chains_dir, 'run?-II-outbreak.rep?.1e6.chains.txt'))
        seqs_fasta = os.path.join(chains_dir, 'II-outbreak-v3.aligned.pruned.fasta')
        seqs_dates = os.path.join(chains_dir, 'II-outbreak-v3.dates.txt')
        print(chains_files)
        pair2info = gather_infector_stats(seqs_fasta, seqs_dates, chains_files, 'infection_pairs.tsv', 
                                          n_burnin=101, pair_reporting_threshold=.1)


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
    #compute_kmc_union()
    #find_linking_kmers_in_kmc_dbs()
    #gather_linkingkmer_stats()
    #make_samples_blastdbs()
    #remove_mumps_contigs()
    #plot_contact_tracing_info()
    compute_num_diffs()

