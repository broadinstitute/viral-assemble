'''
    SPAdes - St. Petersburg Assembler

'''

import logging
import os
import os.path
import subprocess
import shutil
import random

import tools
import tools.samtools
import tools.picard
import util.file
import util.misc

TOOL_NAME = 'spades'
TOOL_VERSION = '3.10.1'

log = logging.getLogger(__name__)

class SpadesTool(tools.Tool):
    """Tool wrapper for SPAdes tool (St. Petersburg Assembler)"""

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, executable='spades.py',
                                                  verifycmd='echo spades.py --test', env='spades_env')]
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

    def assemble(self, reads_fwd, reads_bwd, contigs_out, reads_unpaired = None, contigs_in = None ):
        """Extend contigs using the reads"""

        spades_dir = os.path.join( util.file.make_tmp_dir('spades-dir'), 'spades-out' )
        if 'VNGS_SPADES_DIR' in os.environ: 
            spades_dir = os.environ[ 'VNGS_SPADES_DIR' ] + '-' + str(random.randrange(100000))
        log.info('spades_dir=' + spades_dir)
        args = ['-1', reads_fwd, '-2', reads_bwd ]
        if reads_unpaired: args += [ '-s', reads_unpaired ]
        if contigs_in: args += [ '--trusted-contigs', contigs_in ]
        args += [ '--rna', '-o', spades_dir ]
        
        self.execute( args = args )

        shutil.copyfile( src = os.path.join( spades_dir, 'transcripts.fasta' ), dst = contigs_out)


                      
                      

    

