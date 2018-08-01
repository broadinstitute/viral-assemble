'''
    Tool wrapper for the abyss assembler, including its Sealer gap-closing tool.
'''

import itertools
import functools
import operator
import logging
import os
import os.path
import subprocess
import shlex
import shutil
import tempfile
import time

import Bio.Seq

import tools
import tools.samtools
import util.file
import util.misc

TOOL_NAME = 'abyss'
TOOL_VERSION = '2.1.0'

log = logging.getLogger(__name__)

class AbyssTool(tools.Tool):
    """Tool wrapper for the abyss assembler, including abyss-sealer gap closer."""

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, executable='abyss-sealer')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + list(args)
        log.debug('running abyss tool: ' + ' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def _run_sealer(self, reads, filled, *args, **kwargs):
        # file_args = ('-scaffolds', abspath(scaffolds), '-filled', abspath(filled), 
        #              '-reads', ','.join(map(abspath,reads)))
        more_args = ['--'+arg.replace('_','-')+'='+str(val) for arg, val in kwargs.items()]
        with util.file.tmp_dir('_sealer_run_dir') as sealer_run_dir:
            self.execute(list(args)+more_args+['--fix-errors', '--verbose', '--output-prefix='+os.path.join(sealer_run_dir, 'result')]
                         +list(reads))
            shutil.copyfile(os.path.join(sealer_run_dir, 'result_scaffold.fa'), filled)

    def gapfill(self, in_scaffold, in_bam, out_scaffold, solid_kmer_thresholds=(3,), kmer_sizes=(31, 18, 16),
                min_gap_to_close=4, sealer_opts='', mem_limit_gb=4.0, threads=None, time_soft_limit_minutes=60.0, random_seed=0):
        """Try to fill the gaps in the given scaffold, using the reads.

        Inputs:
            in_scaffold: a FASTA file containing the scaffold.  Each FASTA record corresponds to one
                segment (for multi-segment genomes).  Contigs within each segment are
                separated by Ns.  The exact number of Ns between contigs does not matter, as the length of the gap is one 
                of the things determined by the gap-filling tool.  (But see `min_gap_to_close`).
            in_bam: reads to use for filling the gaps.  Only paired-end reads from the bam file are used, any unpaired reads
                are ignored.
           
        Outputs:
            out_scaffold: the input scaffold, with some of the gaps between contigs possibly filled.

        Params:
            solid_kmer_thresholds: kmers must appear at least this many times in the reads to be considered solid.
                We try gapfilling for all combinations of values of solid_kmer_thresholds and kmer_sizes.
            kmer_sizes: kmer sizes to use.  We try gapfilling for all combinations of values of solid_kmer_thresholds and kmer_sizes.
            min_gap_to_close: stop gap-closing if all gaps are no longer than this many Ns
            sealer_opts: extra command-line flags to pass to Sealer
            mem_limit_gb: max memory to use, in gigabytes
            threads: number of threads to use; None means use all available cores.
            time_soft_limit_minutes: stop trying to close more gaps after this many minutes (currently this is a soft/advisory limit)
            random_seed: random seed for choosing random paths (0 to use current time)
        
        """
        solid_kmer_thresholds = sorted(util.misc.make_seq(solid_kmer_thresholds), reverse=True)
        kmer_sizes = sorted(util.misc.make_seq(kmer_sizes), reverse=True)
        stop_time = time.time() + 60*time_soft_limit_minutes
        threads = util.misc.sanitize_thread_count(threads, tool_max_cores_value=0)
        with tools.samtools.SamtoolsTool().bam2fq_tmp(in_bam) as reads, util.file.tmp_dir('_sealer_dir') as sealer_dir:

            # We call Sealer for a range of parameter combinations.  Output of each call is input to the next call, so
            # each call only deals with gaps not closed by prior calls.  We first try to close using higher-quality kmers,
            # and if that fails try with lower-quality ones.
            prev_scaffold = in_scaffold
            for solid_kmer_threshold, kmer_size in itertools.product(solid_kmer_thresholds, kmer_sizes):

                # here could use kmc to hard-mask from the reads kmers that are not solid

                if not any('N'*min_gap_to_close in str(rec.seq) for rec in Bio.SeqIO.parse(prev_scaffold, 'fasta')):
                    log.info('no gaps left, quittting Sealer early')
                    break
                if time.time() > stop_time:
                    log.info('Time limit for gap closing reached')
                    break

                filled_scaffold = os.path.join(sealer_dir, 'sealer-filled.s{}.k{}.fasta'.format(solid_kmer_threshold, kmer_size))
                self._run_sealer(reads, filled_scaffold,
                                 *shlex.split(sealer_opts),
                                 kmer=kmer_size, bloom_size='5g', input_scaffold=prev_scaffold,
                                 threads=threads)

                prev_scaffold = filled_scaffold

            shutil.copyfile(src=prev_scaffold, dst=out_scaffold)
