'''
    Gap2Seq - assembly gap closing tool
'''

import logging
import os
import os.path
import subprocess
import shutil
import tempfile

import tools
import tools.samtools
import util.file

TOOL_NAME = 'gap2seq'
TOOL_VERSION = '2.0'

log = logging.getLogger(__name__)

class Gap2SeqTool(tools.Tool):
    """Tool wrapper for ABySS assembler."""

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, executable='Gap2Seq.sh')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args, stdout=None):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + args
        log.debug(' '.join(tool_cmd))
        if stdout:
            stdout = open(stdout, 'w')
        subprocess.check_call(tool_cmd, stdout=stdout)
        if stdout:
            stdout.close()

    def gapfill(self, in_scaffold, inBam, out_scaffold, gap2seq_params='', threads=1 ):
        """Try to fill the gaps in the given scaffold, using the reads.

        Args:

           in_scaffold: a FASTA file containing the scaffold.  Contigs within each chromosome are
              separated by Ns.  The number of Ns separating two contigs does not matter.

        """
        samtools = tools.samtools.SamtoolsTool()
        with tempfile.TemporaryDirectory('_gap2seq') as gap2seq_dir:
            with samtools.bam2fq_tmp(inBam) as (reads1, reads2):
                prev_scaffold=in_scaffold
                done=False
                for s in (3,2):
                    for k in (90, 80, 70, 60, 50, 40, 31):
                        gap2seq_filled = os.path.join( gap2seq_dir, 'gap2seq-filled.s{}.k{}.fasta'.format(s,k) )
                        log.info('s={} k={} gap2seq_filled={}'.format(s,k,gap2seq_filled))

                        args = gap2seq_params.split() + [ '-scaffolds', prev_scaffold, '-filled', gap2seq_filled, 
                                                         '-reads', ','.join((reads1,reads2)), '-solid', str(s), 
                                                          '-k', str(k), '-all-upper',
                                                          '-verbose', '-nb-cores', 
                                                         str(threads) ]

                        self.execute( args = args )
                        prev_scaffold=gap2seq_filled
                        prev_scaffold_txt=util.file.slurp_file(prev_scaffold)
                        log.info('gap2seq: s={} k={} Ns={}'.format(s,k,prev_scaffold_txt.count('N')))
                        if 'NNNNNNNNNN' not in prev_scaffold_txt: 
                            log.info('no runs of N left, quittting gap2seq early')
                            done=True
                            break
                    if done: break

                shutil.copyfile( src = prev_scaffold, dst = out_scaffold )

