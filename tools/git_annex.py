'''
    Tool wrapper for git-annex
'''

# * imports

import logging
import collections
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

TOOL_NAME = 'git-annex'
TOOL_VERSION = '7.20181105'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)

# * class GitAnnexTool
class GitAnnexTool(tools.Tool):

    '''Tool wrapper for git-annex.

    TODO:
       - support unlocked annexed files
    '''

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, channel='conda-forge',
                                                  executable='git-annex',
                                                  # install git-annex into own environment, because it has
                                                  # some dependency requirements (like specific Perl version)
                                                  # that conflict with other pkgs.
                                                  env='vngs-git-annex-env')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args, tool='git-annex'):    # pylint: disable=W0221
        """Run program `tool` with arguments `args`."""
        bin_dir = os.path.dirname(self.install_and_get_path())
        tool_cmd = [os.path.join(bin_dir, tool)] + list(map(str, args))

        # ensure that the git and git-annex binaries we installed are first in PATH.
        env = dict(os.environ)
        env['PATH'] = bin_dir + ':' + env['PATH']
        _log.debug(' '.join(tool_cmd))
        _log.debug('PATH=%s', env['PATH'])
        subprocess.check_call(tool_cmd, env=env)

    def execute_git(self, args):
        """Run a git command"""
        self.execute(args, tool='git')

    def print_version(self):
        """Print git-annex version"""
        self.execute(['version'])

    def init_repo(self):
        """Initialize git-annex repo in current dir"""
        self.execute_git(['init'])
        self.execute(['init'])

    def add(self, fname):
        """Add a file to git-annex"""
        self.execute(['add', fname])

    def commit(self, msg):
        """Execute git commit"""
        self.execute_git(['commit', '-m', '"{}"'.format(msg)])

    def initremote(self, name, remote_type, **kw):
        """Initialize a git-annex special remote"""
        remote_attrs = dict(kw)
        if 'encryption' not in remote_attrs: remote_attrs['encryption'] = 'none'
        self.execute(['initremote', name, 'type='+remote_type]+['='.join((k,v)) for k, v in remote_attrs.items()])

    def move(self, fname, to_remote_name):
        """Move a file from local repo to a remote"""
        self.execute(['move', fname, '--to', to_remote_name])

    def _get_link_into_annex(self, f):
        """If `f` points to an annexed file, possibly through a chain of symlinks, return
        information about the symlink actually pointing into the annex (details below).
        If `f` does not point to an annexed file, return (f, None).

        Details: if `f` points into the annex through a chain of symlinks, returns the pair
        (link_into_annex, target_of_link_into_annex), where link_into_annex is the path to the
        symlink that points directly into the annex (i.e. the penultimate link in a chain of symlinks),
        while target_of_link_into_annex is the path within the annex to which link_into_annex points.

        TODO: handle unlocked git-annex files.
        """
        annex_link_target = None
        f_cur = f
        paths_seen = set()
        while os.path.islink(f_cur):
            _log.debug('_get_link_into_annex: f_cur=%s cwd=%s', f_cur, os.getcwd())
            f_cur_abs = os.path.abspath(f_cur)
            if f_cur_abs in paths_seen:
                _log.warning('Circular symlinks! %s', f_cur_abs)
                break
            paths_seen.add(f_cur_abs)

            link_target = os.readlink(f_cur)
            if '.git/annex/objects/' in link_target:
                annex_link_target = link_target
                f = f_cur
                break
            if os.path.isabs(link_target):
                f_cur = link_target
            else:
                f_cur = os.path.join(os.path.dirname(f_cur), link_target)
        return f, annex_link_target

    def is_annexed_file(self, f):
        """Test is `f` is a file controlled by git-annex."""
        return os.path.lexists(f) and not os.path.isdir(f) and self._get_link_into_annex(f)[1]

    def lookupkey(self, f):
        """Get the git-annex key of an annexed file.  Note that, unlike git-annex-lookupkey, this looks at the file in the
        working copy, not in the index."""
        
        link_into_annex, target_of_link_into_annex = self._get_link_into_annex(f)
        util.misc.chk(target_of_link_into_annex)
        return os.path.basename(target_of_link_into_annex)

    def examinekey(self, key):
        """Return a dict of info that can be gleaned from the key itself"""
        pass

    def get_annexed_file_attrs(self, f):
        """Get annexed file info that can be determined quickly e.g. from the key"""
        raise NotImplemented()


    def get(self, f):
        """Ensure the file exists in the local annex, fetching it from a remote if necessary.
        Unlike git-annex-get, follows symlinks and  will get the file regardless of what the current dir is."""

        # TODO: what if f is a dir, or list of dirs?  including symlinks?
        # then, to preserve git-annex-get semantics, need to walk through the files under the tree.

        _log.debug('git-annex get %s (cwd=%s)', f, os.getcwd())
        assert os.path.islink(f)
        f, link_target = self._get_link_into_annex(f)
        if not os.path.isfile(f):
            with util.file.pushd_popd(os.path.dirname(os.path.abspath(f))):
                self.execute(['get', os.path.basename(f)])
        assert os.path.isfile(f)

    def drop(self, f):
        """Drop the file from its local repo."""
        _log.debug('git-annex drop %s (cwd=%s)', f, os.getcwd())
        assert os.path.islink(f)
        f, link_target = self._get_link_into_annex(f)
        if os.path.isfile(f):
            with util.file.pushd_popd(os.path.dirname(os.path.abspath(f))):
                self.execute(['drop', os.path.basename(f)])
        assert os.path.islink(f)
        assert not os.path.isfile(f)

# end: class GitAnnexTool
