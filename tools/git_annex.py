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
import functools
import contextlib
import copy

import tools
import util.file
import util.misc

TOOL_NAME = 'git-annex'
TOOL_VERSION = '7.20181211'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)


# class _Batcher(object):
#     """Batches """

# def _batcher(single_method, batch_method):
#     single_method.argpacks = []
#     def single_method_impl(self, *args):
#         if not self._batching:
#              batch_method([args])
#         single_method.argpacks.append(args)

# * class GitAnnexTool
class GitAnnexTool(tools.Tool):

    '''Tool wrapper for git-annex.

    '''

    # Implementation notes

    # Fields:
    #   _batching: if True, commands are batched rather then run immediately
    #   _batched_cmds: map from git-annex command to input tuples to feed to the command in batch mode

    # TODO:
    #    - support unlocked annexed files
    #    - add waiting for git index.lock to clear for some ops?  which ones?
    #    - automatically decide what to batch together?
    #    - support caching of metadata?
    #    - check that batched commands are indeed independent?  (put names of affected repo paths in shared mem?)
    #    - (option to) do opos on a separate checkout in a tmp  dir and then merge?  checkout first commit if only adding dirs

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, channel='conda-forge')]
        tools.Tool.__init__(self, install_methods=install_methods)
        self._batching = False
        self._batched_cmds = collections.defaultdict(list)

    def version(self):
        return TOOL_VERSION

    def _get_bin_dir(self):
        return os.path.dirname(self.install_and_get_path())

    def _get_run_env(self):
        if not hasattr(self, '_run_env'):
            # ensure that the git and git-annex binaries we installed are first in PATH.
            self._run_env = dict(os.environ, PATH=':'.join((self._get_bin_dir(), os.environ['PATH'])))
        return self._run_env

    def _call_tool_cmd(self, tool_cmd, **kw):
        subprocess.check_call(tool_cmd, env=self._get_run_env(), **kw)

    def execute(self, args, tool='git-annex', batch_args=()):    # pylint: disable=W0221
        """Run program `tool` with arguments `args`."""

        tool_cmd = (os.path.join(self._get_bin_dir(), tool),) + tuple(map(str, args))

        _log.debug(' '.join(tool_cmd))

        if not batch_args:
            _log.debug('Instant non-batched (no batch args): %s', tool_cmd)
            self._call_tool_cmd(tool_cmd)
            return

        self._batched_cmds[tool_cmd].append(batch_args)
        if not self._batching:
            _log.debug('Instantly running cmd %s with batch args %s', tool_cmd, batch_args)
            self.execute_batched_commands(instant=True)

    def execute_batched_commands(self, instant=False):
        """Run any saved batched commands.

        Note that if one command fails, others won't be run.
        """

        # TODO: figure out correct behavior if some cmds fail

        _log.debug('RUNNING BATCH CMDS: instant=%s', instant)
        while self._batched_cmds:
            tool_cmd, batch_inputs = self._batched_cmds.popitem()
            _log.debug('Running batch cmd %s with batched args %s', tool_cmd, batch_inputs)
            with util.file.tempfname(prefix='git-annex-batch-inputs') as batch_inps_fname:
                util.file.dump_file(batch_inps_fname,
                                    '\n'.join([' '.join(map(str, inps)) for inps in batch_inputs]))
                with open(batch_inps_fname) as batch_inps_file:
                    _log.debug('CALLING BATCH: %s batch_inps_file=%s', tool_cmd, batch_inps_file)
                    subprocess.check_call(tool_cmd, stdin=batch_inps_file)

    @contextlib.contextmanager
    def batching(self):
        """Commands run within this context will, on successful context exit, be run as a batch (if not already batching) or
        added to currently accumulating batch (if already batching).  If an exception occurs within the context,
        commands batched within the context will not be run or saved for later.
        """

        # TODO: figure out correct behavior if exception occurs within context
        
        _log.debug('BATCHING enter context: self._batching was %s', self._batching)
        save_batching = self._batching
        save_batched_cmds = copy.deepcopy(self._batched_cmds)
        self._batching = True

        try:
            yield
            _log.debug('BATCHING: successfully returned from yield')
            self._batching = save_batching
            if not self._batching:
                self.execute_batched_commands()
        except Exception as e:
            _log.debug('BATCHING: exception %s', e)
            self._batching = save_batching
            self._batched_cmds = save_batched_cmds
            raise

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
        _log.debug('CALL TO ADD %s; batch status = %s', fname, self._batching)
        self.execute(['add', '--batch'], batch_args=(fname,))
        _log.debug('RETURNED FROM CALL TO ADD %s; batch status = %s', fname, self._batching)

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

    def fromkey(self, key, fname):
        """Manually set up a symlink to a given key as a given file"""
        self.execute(['fromkey', '--force', '--batch'], batch_args=(key, fname))

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

    def is_link_into_annex(self, f):
        """Test if `f` is a file controlled by git-annex (not necessarily present in local repo)."""
        return os.path.lexists(f) and not os.path.isdir(f) and self._get_link_into_annex(f)[1]

    def is_file_in_annex(self, f):
        """Tests if `f` (possibly at the end of a chain of symlinks) is a file present in the local annex."""
        return self.is_link_into_annex(f) and os.path.isfile(f)

    def lookupkey(self, f):
        """Get the git-annex key of an annexed file.  Note that, unlike git-annex-lookupkey, this looks at the file in the
        working copy, not in the index."""
        
        link_into_annex, target_of_link_into_annex = self._get_link_into_annex(f)
        util.misc.chk(target_of_link_into_annex)
        return os.path.basename(target_of_link_into_annex)

    def examinekey(self, key):
        """Return a dict of info that can be gleaned from the key itself.
        Most accurate would be to run 'git annex examinekey', but for speed we do the
        parsing ourselves.
        """
        attrs = collections.OrderedDict()
        key_attrs_str, key_name = key.rsplit('--', 1)
        attrs['key_name'] = key_name
        if key.startswith('MD5'):
            attrs['md5'] = key_name[:32]
            key_attrs_parts = key_attrs_str.split('-')
            key_attrs_size_part = [p for p in key_attrs_parts if p.startswith('s')][0]
            attrs['size'] = int(key_attrs_size_part[1:])
        return attrs

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
