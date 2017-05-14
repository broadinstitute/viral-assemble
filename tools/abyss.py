'''
    ABySS - the ABySS genome assembler, including the Sealer tool for filling gaps between contigs.
'''

import logging
import os
import os.path
import subprocess
import shutil
import tempfile

import tools
import tools.samtools

TOOL_NAME = 'abyss'
TOOL_VERSION = '2.0.1'

log = logging.getLogger(__name__)

class AbyssTool(tools.Tool):
    """Tool wrapper for ABySS assembler."""

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, executable='abyss-sealer')]
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

    def gapfill(self, in_scaffold, inBam, out_scaffold, sealer_params, threads=1 ):
        """Try to fill the gaps in the given scaffold, using the reads.

        Args:

           in_scaffold: a FASTA file containing the scaffold.  Contigs within each chromosome are
              separated by Ns.  The number of Ns separating two contigs does not matter.

        """
        samtools = tools.samtools.SamtoolsTool()
        with tempfile.TemporaryDirectory('_sealer') as sealer_dir:
            with samtools.bam2fq_tmp(inBam) as (reads1, reads2):
                sealer_pfx = os.path.join( sealer_dir, 'sealer-out' )

                args = sealer_params.split() + [ '-S', in_scaffold, '-o', sealer_pfx, '-j', str(threads), '-v',
                                                 reads1, reads2 ]

                self.execute( args = args )

                shutil.copyfile( src = sealer_pfx + '_scaffold.fa', dst = out_scaffold )
                #shutil.copyfile( src = sealer_pfx + '_merged.fa', dst = out_scaffold + '.merged.fa' )
            


                      
                      

    

