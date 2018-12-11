'''
    Tool wrapper for the Cromwell workflow manager.
'''

import logging
import os
import os.path
import subprocess
import shutil
import random
import shlex
import tempfile
import pipes
import time

import tools
import util.file
import util.misc

TOOL_NAME = 'cromwell'
TOOL_VERSION = '0.36'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)

class CromwellTool(tools.Tool):

    '''Tool wrapper for Cromwell workflow manager'''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION)]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + list(map(str, args))
        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    

