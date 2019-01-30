'''
    Tool wrapper for the DNAnexus CLI.
'''

# * imports

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
import platform

import tools
import util.file
import util.misc

TOOL_NAME = 'dxpy'
TOOL_VERSION = '0.273.0'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)

# * class DxTool

class DxTool(tools.Tool):

    '''Tool wrapper for DNAnexus CLI'''

# ** init, execute
    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION,
                                                  executable='dx')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + list(map(str, args))
        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def print_version(self):
        """Print dx cli version"""
        self.execute(['--version'])

# * end
# end: class DxTool(tools.Tool)

