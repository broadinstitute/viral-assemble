'''
    Tool wrapper for docker.
'''

__all__ = ['DockerTool']

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
import re

import tools
import util.file
import util.misc

TOOL_NAME = 'docker'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)

# * Misc utils
def _noquote(s):
    return '_noquote:' + str(s)

def _quote(s):
    s = str(s)
    if s.startswith('_noquote:'): return s[len('_noquote:'):]
    return pipes.quote(s) if hasattr(pipes, 'quote') else shlex.quote(s)

def _make_cmd(cmd, *args):
    _log.debug('ARGS=%s', args)
    return ' '.join([cmd] + [_quote(str(arg)) for arg in args if arg not in (None, '')])

def _run(cmd, *args):
    cmd = _make_cmd(cmd, *args)
    _log.debug('running command: %s cwd=%s', cmd, os.getcwd())
    beg_time = time.time()
    subprocess.check_call(cmd, shell=True)
    _log.debug('command succeeded in {}s: {}'.format(time.time()-beg_time, cmd))

def _run_get_output(cmd, *args):
    cmd = _make_cmd(cmd, *args)
    _log.debug('running command: %s cwd=%s', cmd, os.getcwd())
    beg_time = time.time()
    output = subprocess.check_output(cmd, shell=True)
    _log.debug('command succeeded in {}s: {}'.format(time.time()-beg_time, cmd))
    return output.strip()

# * class DockerTool

class DockerTool(tools.Tool):

    '''Tool wrapper for Docker'''

# ** init, execute
    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.PrexistingUnixCommand(path=distutils.spawn.find_executable('docker'))]
        tools.Tool.__init__(self, install_methods=install_methods)

    @staticmethod
    def is_docker_installed():
        return distutils.spawn.find_executable('docker')

    def version(self):
        return TOOL_VERSION

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + list(map(str, args))
        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def has_image_hash(self, docker_img):
        """Test whether the given docker_img includes the sha256 hash"""
        return re.search(r'@sha256:[0-9a-z]{64}\Z', docker_img) 

    def add_image_hash(self, docker_img):
        """Append to `docker_img` the sha256 hash of the image, if not yet present."""
        return docker_img + ('' if self.has_image_hash(docker_img) else '@'+self.get_docker_hash(docker_img))

    def strip_image_hash(self, docker_img):
        """Strip from `docker_img` the sha256 image hash, if present."""
        if self.has_image_hash(docker_img):
            docker_img = docker_img[:-(len('@sha256:')+64)]
        return docker_img
        
    def get_docker_hash(self, docker_img):
        """Return a docker hash, given a docker tag."""
        _run('docker pull ' + docker_img)
        if docker_img.startswith('sha256:'):
            return docker_img
        digest_lines = _run_get_output('docker', 'images', '--digests', '--no-trunc', '--format',
                                       '{{.Repository}}:{{.Tag}} {{.Digest}}', _noquote('|'), 'grep',
                                       docker_img+(':' if ':' not in docker_img else '') + ' sha256:')
        digest_lines = [line for line in util.misc.maybe_decode(digest_lines).rstrip('\n').split('\n') if docker_img in line]
        assert len(digest_lines) == 1
        digest_line = digest_lines[0]
        _log.debug('digest_line is |||{}|||'.format(digest_line))
        img, digest = digest_line.strip().split()
        assert img == docker_img + (':latest' if ':' not in docker_img else '')
        assert digest.startswith('sha256:') and len(digest) == 71
        return digest

# * end
# end: class DockerTool(tools.Tool)


