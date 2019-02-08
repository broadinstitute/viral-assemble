'''
    Tool wrapper for docker.
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
import distutils.spawn

import tools
import util.file
import util.misc

TOOL_NAME = 'docker'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)

# * class DockerTool

class DockerTool(tools.Tool):

    '''Tool wrapper for Docker'''

# ** init, execute
    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PreexistingUnixCommand(path=distutils.spawn.find_executable('docker'))]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + list(map(str, args))
        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

# * end
# end: class DockerTool(tools.Tool)


