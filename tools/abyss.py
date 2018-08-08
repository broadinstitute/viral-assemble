'''
    ABySS - the ABySS genome assembler.
'''

import logging
import os
import os.path
import subprocess
import shutil
import random
import tempfile

import tools
import tools.samtools
import tools.picard
import util.file
import util.misc

TOOL_NAME = 'abyss'
TOOL_VERSION = '2.1.0'

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

    def gapfill(self, in_scaffold, reads_fwd, reads_bwd, out_scaffold, kmer_sizes, bloom_filter_size_gb,
                max_paths=4,
                fix_errors=False, n_threads=1, sealer_opts=[] ):
        """Extend contigs using the reads"""

        with tempfile.TemporaryDirectory('_sealer') as sealer_dir:
            sealer_pfx = os.path.join( sealer_dir, 'sealer-out' )
            log.debug('sealer_out=' + sealer_pfx)
            args = [ '-v', '-b'+ str(bloom_filter_size_gb) + 'g' ] + [ '-k' + str(k) for k in kmer_sizes ] + \
                   [ '-S', in_scaffold, '-P', str(max_paths)  ]
            args += [ '-o', sealer_pfx ]
            args += sealer_opts
            args += [ reads_fwd, reads_bwd ]

            self.execute( args = args )

            shutil.copyfile( src = sealer_pfx + '_scaffold.fa', dst = out_scaffold )
            shutil.copyfile( src = sealer_pfx + '_merged.fa', dst = out_scaffold + '.merged.fa' )



                      
                      

    

