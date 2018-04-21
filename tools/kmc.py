'''
    Tool wrapper for the KMC kmer counter 
    ( http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about )
'''

import logging
import os
import os.path
import subprocess
import shutil
import shlex

import tools
import util.file
import util.misc

TOOL_NAME = 'kmc'
TOOL_VERSION = '3.0.0'

log = logging.getLogger(__name__)

class KmcTool(tools.Tool):
    '''Tool wrapper for KMC kmer counter'''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, executable='kmc',
                                                  verifycmd='kmc')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def build_kmer_db(self, args):    # pylint: disable=W0221
        """Build a database of kmers."""
        with util.file.tmp_dir(suffix='kmcdb') as t_dir:
            tool_cmd = [self.install_and_get_path()] + list(map(str, args)) + [t_dir]
            log.info(' '.join(tool_cmd))
            print('tool_cmd=', tool_cmd)
            subprocess.check_call(tool_cmd)

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()+'_tools'] + list(map(str, args))
        log.info(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

# end: class KmcTool(tools.Tool)

