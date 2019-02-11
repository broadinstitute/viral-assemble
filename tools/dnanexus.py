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

try:
    from urllib import urlencode, pathname2url
except ImportError:
    from urllib.parse import urlencode
    from urllib.request import pathname2url

import dxpy
import uritools

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

    DX_URI_PFX='dx://'

    @classmethod
    def url_to_dxid(cls, url):
        """Given a dx://file-xxx URL, return the dxid in it"""
        return uritools.urisplit(url).authority
    
    @classmethod
    def dx_make_download_url(cls, dxid, duration='2h'):
        return _run_get_output('dx', 'make_download_url', dxid, '--duration', duration)

    @classmethod
    def standardize_dx_url(cls, url):
        dxid = _url_to_dxid(url)
        dx_descr = _dx_describe(_url_to_dxid(url))
        return cls.DX_URI_PFX + dxid + '/' + pathname2url(dx_descr['name'])

#    @util.misc.memoize_persist(to_picklable=functools.partial(json.dumps, separators=(',',':')),
#                               from_picklable=_json_loads)
    @staticmethod
    def describe(dxid):
        """Return json description for the given dxid"""
        return dxpy.describe(dxid)
        #return _run_get_json('dx', 'describe', '--verbose', '--details', '--json', dxid)

# * end
# end: class DxTool(tools.Tool)
