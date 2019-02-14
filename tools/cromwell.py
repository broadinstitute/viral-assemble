'''
    Tool wrapper for the Cromwell workflow manager.
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
import contextlib

import cromwell_tools
import cromwell_tools.cromwell_api
import cromwell_tools.cromwell_auth

import tools
import util.file
import util.misc

TOOL_NAME = 'cromwell'
TOOL_VERSION = '0.36'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)

# * class CromwellTool

class CromwellTool(tools.Tool):

    '''Tool wrapper for Cromwell workflow manager'''

# ** init, execute
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

    class CromwellServer(object):
        """Represents a specific running cromwell server'"""

        def __init__(self, cromwell_tool, auth):
            self.cromwell_tool = cromwell_tool
            self.auth = auth
            self.api = cromwell_tools.cromwell_api.CromwellAPI()
            self.cromwell_process = subprocess.Popen([cromwell_tool.install_and_get_path(), 'server'])

        def shutdown(self, timeout=300):
            """Shut down the cromwell server"""
            self.cromwell_process.terminate()
            self.cromwell_process.wait(timeout=timeout)

        def health(self, *args, **kwargs):
            """Do nothing is the server is running fine, else raise a RuntimeError"""
            return self.api.health(self.auth, *args, **kwargs)

        def is_healthy(self):
            health_response = self.health()
            return health_response and not (400 <= health_response.status_code < 600)

    @contextlib.contextmanager
    def cromwell_server(self):
        url = 'http://localhost:8000'
        server = self.CromwellServer(cromwell_tool=self,
                                     auth=cromwell_tools.cromwell_auth.CromwellAuth.from_no_authentication(url=url))
        util.misc.chk(server.is_healthy())
        yield server
        server.shutdown()

# ** Metadata handling


# * end
# end: class CromwellTool(tools.Tool)
