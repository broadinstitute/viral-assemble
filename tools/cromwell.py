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
import collections
import time
import signal

import cromwell_tools
import cromwell_tools.cromwell_api
import cromwell_tools.cromwell_auth

import tools
import util.file
import util.misc

TOOL_NAME = 'cromwell'
TOOL_VERSION = '0.37'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)

# * class CromwellTool

class CromwellTool(tools.Tool):

    '''Tool wrapper for Cromwell workflow manager'''

# ** init, execute
    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, env='vngs_cromwell_env')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + list(map(str, args))
        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    class CromwellServer(object):
        """Represents a specific running cromwell server'"""

        def __init__(self, cromwell_tool, port=8000):
            self.cromwell_tool = cromwell_tool
            self.url = 'http://localhost:{}'.format(port)
            self.auth = cromwell_tools.cromwell_auth.CromwellAuth.from_no_authentication(url=self.url)
            self.api = cromwell_tools.cromwell_api.CromwellAPI()
            args = [cromwell_tool.install_and_get_path(), '-Dwebsevice.port={}'.format(port), 'server']
            _log.info('starting cromwell server: args=%s auth=%s', args, self.auth)
            self.cromwell_process = subprocess.Popen(args)
            time.sleep(2)

        def shutdown(self, timeout=300):
            """Shut down the cromwell server"""
            util.misc.kill_proc_tree(self.cromwell_process)
            _log.info('Waiting for Cromwell to terminate, timeout=%d', timeout)
            self.cromwell_process.wait(timeout=timeout)
            _log.info('Cromwell terminated successfully')

        def health(self, *args, **kwargs):
            """Do nothing is the server is running fine, else raise a RuntimeError"""
            return self.api.health(self.auth, *args, **kwargs)

        def is_healthy(self):
            health_response = self.health()
            if health_response.status_code != 200:
                return False
            try:
                health_report = util.misc.json_loads(health_response.content)
            except Exception:
                return False
            return isinstance(health_report, collections.Mapping) and \
                all(status.get('ok', False) for subsystem, status in health_report.items())

    @contextlib.contextmanager
    def cromwell_server(self, port=8000):
        """Start a cromwell server, shut it down when context ends."""
        server = self.CromwellServer(cromwell_tool=self, port=port)
        time.sleep(10)
        _log.info('IN CROMWELL, AUTH IS %s HLTH IS %s', server.auth, server.health())
        util.misc.chk(server.is_healthy())
        try:
            yield server
        finally:
            server.shutdown()

# ** Metadata handling


# * end
# end: class CromwellTool(tools.Tool)
