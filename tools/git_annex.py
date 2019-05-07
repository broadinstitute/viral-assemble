'''
    Tool wrapper for git-annex (and git)
'''

# * imports

from __future__ import print_function
import logging
import collections
import os
import os.path
import sys
import stat
import concurrent.futures
import uuid
import contextlib2
import random

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

import shutil
import random
import shlex
import tempfile
import pipes
import time
import functools
import operator
import copy
import contextlib2 as contextlib
import uritools
import threading
import dxpy

import tools
import tools.gcloud
import tools.dnanexus
import util.file
import util.misc
import util.version

TOOL_NAME = 'git-annex'
TOOL_VERSION = '7.20190507'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)

class _ValHolder(object):

    """A mutable holder for one value.  Used to provide a place to record output of commands."""

    def __init__(self, val=None):
        self.val = val

    def __call__(self, val):
        self.val = val

# NamedTuple: _Batch - a batch of calls to git-annex, involving one git-annex command with a list of different arg combos
#   Fields:
#      priority: a priority value; batches with higher priority will be run earlier
#      tool_cmd: the invariant part of the command
#      cwd: the current directory in which the batch should be invoked
#      preproc: an optional function to call to preprocess the batch.
_Batch = collections.namedtuple('_Batch', ['priority', 'tool_cmd', 'cwd', 'preproc'])

# NamedTuple: _BatchedCall - arguments to one git-annex call in a set of batched cals, and an optional place to record call output.
_BatchedCall = collections.namedtuple('_BatchedCall', ['batch_args', 'output_acceptor'])



# * class GitAnnexTool
class GitAnnexTool(tools.Tool):

# ** class docs

    '''Tool wrapper for git-annex.

    Notes:

    The repository on which the commands operate is given implicitly by the current working directory.
    (This may change later -- passing the repository explicitly an as object may be a better design.)


    There is support for batching of operations.  You can start a batching context by calling the batching() method,
    which will return a "batching" version of this GitAnnexTool.  When you invoke a command on this batching version,
    you can indicate that it need not run immediately, by adding now=False to the argument list; the command then may be
    run anytime between immediately and the end of the batching context.  

    '''

    # Implementation notes

    # Fields:
    #   _batched_cmds: None or map from git-annex command to list of _BatchedCall records

    # TODO:
    #    - pass -J options or configure annex.jobs
    #    - support unlocked annexed files
    #    - add waiting for git index.lock to clear for some ops?  which ones?
    #    - automatically decide what to batch together?
    #    - support caching of metadata?
    #    - check that batched commands are indeed independent?  (put names of affected repo paths in shared mem?)
    #    - (option to) do opos on a separate checkout in a tmp  dir and then merge?  checkout first commit if only adding dirs
    #    - use pyfilesystem2 to deal with the different filesystems, including staging.
    #    - decorator to do changes on a temp branch and only merge if succeed; branch the git_annex branch as well.

# ** general infrastructure for running commands

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, channel='conda-forge')]
        tools.Tool.__init__(self, install_methods=install_methods)
        self._batched_cmds = None
        self._url2filestat = collections.OrderedDict()
        self._executor = None
        _log.debug('Created GitAnnexTool %s', id(self))

    def version(self):
        return TOOL_VERSION

    def _get_bin_dir(self):
        """Return the directory containing the git-annex binary"""
        if not hasattr(self, 'save_bin_dir'):
            bin_dir = os.path.dirname(self.install_and_get_path())
            util.misc.chk(self.is_installed())
            self.save_bin_dir = bin_dir
        return self.save_bin_dir

    def _get_run_env(self):
        """Return the environment in which git-annex should be run.  This environment will include in its PATH
        the binaries for git-annex external special remote implementations that we add."""
        if not hasattr(self, '_run_env'):
            # ensure that the git and git-annex binaries we installed are first in PATH.
            self._run_env = dict(os.environ,
                                 PATH=':'.join((self._get_bin_dir(),
                                                os.path.join(util.version.get_project_path(),
                                                             'tools', 'git-annex-remotes'),
                                                os.environ['PATH'])),
                                 PYTHONPATH=':'.join((util.version.get_project_path(), os.environ.get('PYTHONPATH', ''))))
        return self._run_env

    @contextlib.contextmanager
    def maybe_now(self, now=True):
        """Start a new batching context if `now` is True, else do nothing to inherit an existing batching context."""
        if now:
            with self.batching() as self:
                yield self
        else:
            yield self

    def _add_now_arg(method):
        """To `method`, add a keyword-only argument `now`, such that when it is True (defaut), the method runs immediately,
        while if it is False, the method runs when the current batching context ends."""
        @util.misc.wraps(method)
        def impl(self, *args, **kw):
            with self.maybe_now(kw.get('now', True)) as self:
                kw = {k: v for k, v in kw.items() if k != 'now'}
                return method(self, *args, **kw)
        return impl

    def execute(self, args, batch_args=None, output_acceptor=None, now=True, tool='git-annex', priority=0, preproc=None):
        """Run the given git or git-annex command.  If batching is in effect (see self.batching()), the command is stored and will be
        executed when the batching context exits.
        
        Args:
            args: list of arguments too the command
            batch_args: if given, tuple of args to a batch git-annex command
            output_acceptor: if given, the output of the command will be passed to this callable.
              Note that if batching is in effect, this will only happen when the batching context closes.
            now: execute instantly even if batching is in effect
            tool: git-annex or git, defaults to git-annex
            priority: when executing batched commands, commands with higher priority will be executed earlier
            preproc: if given, called to preprocess a batch
        """
        _log.debug('EXECUTE: args=%s now=%s batch_args=%s self.is_batching=%s', args, now, len(batch_args or ()),
                   self.is_batching())
        for batch_arg in (batch_args or ()):
            _log.debug('batch_arg: %s', batch_arg)
        with self.maybe_now(now=now or not batch_args or not self.is_batching()) as self:
            self._add_command_to_batch(args=args, batch_args=batch_args, output_acceptor=output_acceptor, tool=tool,
                                       priority=priority, preproc=preproc)

    def execute_batch(self, args, batch_args, output_acceptor=None, priority=0, preproc=None):
        """Record a git-annex command to be run at the end of the current batching context."""
        return self.execute(args=args, batch_args=batch_args, output_acceptor=output_acceptor, now=False, priority=priority,
                            preproc=preproc)

    def is_batching(self):
        """Return True if batching is in effect"""
        return self._batched_cmds is not None

    def _add_command_to_batch(self, args, batch_args, output_acceptor, tool, priority, preproc=None):
        """Add command to current batch."""
        if batch_args:
            util.misc.chk(tool == 'git-annex')
            args = tuple(args) + ('--batch',)
        tool_cmd = (os.path.join(self._get_bin_dir(), tool),) + tuple(map(str, args))
        batch = _Batch(priority=priority, tool_cmd=tool_cmd, cwd=os.getcwd(), preproc=preproc)
        batched_call = _BatchedCall(batch_args=batch_args, output_acceptor=output_acceptor)
        self._batched_cmds.setdefault(batch, []).append(batched_call)

    def _execute_batched_commands(self):
        """Run any saved batched commands.

        Note that if one command fails, others won't be run.
        """

        # TODO: figure out correct behavior if some cmds fail

        _log.debug('RUNNING BATCH CMDS: %s', self._batched_cmds)
        for batch in sorted(self._batched_cmds, key=operator.attrgetter('priority'), reverse=True):
            batch_calls = self._batched_cmds.pop(batch)
            # if any are delayed calls, call them now
            if batch_calls:
                def _resolve_futures(args):
                    return [util.misc.maybe_wait_for_result(arg, timeout=1500) for arg in args]
                batch_calls = [batch_call._replace(batch_args=_resolve_futures(batch_call.batch_args)) \
                               if batch_call.batch_args else batch_call \
                               for batch_call in batch_calls]

            if batch.preproc:
                batch_calls = batch.preproc(batch_calls)
            if not batch_calls: continue
            subprocess_call_args = {}
            if batch_calls[0].batch_args:
                # This is a git-annex command with the --batch flag: prepare an input string of the
                # command arguments, one line per batched call.
                _log.debug('calling batch calls: %s %s', batch.tool_cmd, batch_calls)
                

                batch_inp_str = '\n'.join([' '.join(map(str, batch_call.batch_args)) for batch_call in batch_calls]) + '\n'
                subprocess_call_args.update(input=batch_inp_str)

            if batch_calls[0].output_acceptor:
                subprocess_call_args.update(stdout=subprocess.PIPE)

            _log.debug('CALLING SUBPROCESS RUN: %s %s', batch.tool_cmd, subprocess_call_args)
            try:
                result = subprocess.run(batch.tool_cmd, check=True, cwd=batch.cwd, universal_newlines=True,
                                        env=self._get_run_env(), **subprocess_call_args)
                _log.debug('RETURNED FROM SUBPROCESS RUN: %s %s %s', batch.tool_cmd, subprocess_call_args, result)
            except Exception as e:
                _log.debug('FAILED SUBPROCESS RUN: %s %s %s', e, batch.tool_cmd, subprocess_call_args)
                raise

            if batch_calls[0].output_acceptor:
                output = util.misc.maybe_decode(result.stdout).rstrip('\n')
                call_outputs = output.split('\n') if batch_calls[0].batch_args else [output]
                util.misc.chk(len(call_outputs) == len(batch_calls))
                for batch_call, call_output in zip(batch_calls, call_outputs):
                    batch_call.output_acceptor(call_output)
        # end: for batch in sorted(self._batched_cmds, key=operator.attrgetter('priority'), reverse=True):
        util.misc.chk(not self._batched_cmds)
    # end: def _execute_batched_commands(self):

    @contextlib.contextmanager
    def batching(self):
        """Create a batching context.  Commands issued within the batching context _may_ be delayed and issued in batches
        for efficiency; we guarantee that when the batching context exits, the commands have been executed.

        The value returned by this context manager is a GitAnnexTool object that batches commands issued to it,
        and will run them by the time of the context exit.
        """
        _log.debug('ENTERING BATCHING CONTEXT; was batching? {}'.format(self.is_batching()))
        self = copy.copy(self)
        self._batched_cmds = collections.OrderedDict()
        with concurrent.futures.ThreadPoolExecutor(max_workers=util.misc.available_cpu_count()) as self._executor:
            yield self
            _log.debug('EXITING BATCHING CONTEXT; cmds are: {}'.format(self._batched_cmds))
            self._execute_batched_commands()

    def execute_git(self, args, **kw):
        """Run a git command"""
        self.execute(args, tool='git', **kw)

    def execute_and_get_output(self, args, **kw):
        """Run a git or git-annex command and return its output"""
        output_acceptor = _ValHolder()
        self.execute(args, output_acceptor=output_acceptor, **kw)
        return util.misc.maybe_decode(output_acceptor.val).rstrip('\n')

    def execute_git_and_get_output(self, args, **kw):
        """Run a git or git-annex command and return its output"""
        return self.execute_and_get_output(args, tool='git', **kw)

    def get_git_config(self):
        """Return git configuration, as a dict from config param to value.  Note that value is always a string."""
        output_lines = self.execute_git_and_get_output(['config', '-l']).split('\n')
        return collections.OrderedDict(line.rsplit('=', 1) for line in output_lines)

    def get_repo_root(self):
        """Get the root of the repository corresponding to the current dir"""
        return self.execute_git_and_get_output(['rev-parse', '--show-toplevel'])

# ** specific commands

    def print_version(self):
        """Print git-annex version"""
        self.execute(['version'])

    def init_repo(self):
        """Initialize git-annex repo in current dir"""
        self.execute_git(['init'])
        self.execute(['init'])

    def initremote_external(self, remote_name, externaltype):
        """Init an external special remote ldir"""
        self.execute(['initremote', remote_name, 'type=external', 'externaltype={}'.format(externaltype), 'encryption=none'])

# *** git commands

    def commit(self, msg):
        """Execute git commit"""
        self.execute_git(['commit', '-m', '"{}"'.format(msg)])

    def get_first_commit(self):
        """Return sha1 hash of first commit in repo (error if more than one)"""
        parentless_commits = self.execute_git_and_get_output(['rev-list', '--max-parents=0', 'HEAD']).strip().split('\n')
        util.misc.chk(len(parentless_commits) == 1)
        return parentless_commits[0]

    @contextlib.contextmanager
    def tmp_worktree(self, directory, branch, start_branch):
        """Make a temp worktree in a given dir"""
        directory = os.path.abspath(directory)
        with util.file.tmp_dir(dir=directory, suffix='_worktree') as temp_worktree_dir:
            self.execute_git(['worktree', 'add', '-b', branch, temp_worktree_dir, start_branch])
            try:
                with util.file.pushd_popd(temp_worktree_dir):
                    yield temp_worktree_dir
            finally:
                if not util.file.keep_tmp():
                    try:
                        self.execute_git(['worktree', 'remove', temp_worktree_dir])
                    except Exception:
                        _log.warning('Could not remove temp worktree %s', temp_worktree_dir)

    @contextlib.contextmanager
    def _tmp_worktree(self, directory, branch, start_branch):
        """Make a temp worktree in a given dir"""
        directory = os.path.abspath(directory)
        with util.file.tmp_dir(dir=directory, suffix='_worktree') as temp_worktree_dir:
            self.execute_git(['worktree', 'add', '-b', branch, temp_worktree_dir, start_branch])
            try:
                with util.file.pushd_popd(temp_worktree_dir):
                    yield temp_worktree_dir
                    
            finally:
                if not util.file.keep_tmp():
                    try:
                        self.execute_git(['worktree', 'remove', temp_worktree_dir])
                    except Exception:
                        _log.warning('Could not remove temp worktree %s', temp_worktree_dir)

    @contextlib.contextmanager
    def _in_tmp_worktree(self, worktree_group_id='', chdir=True, keep=False):
        """Create a temp branch in a temp worktree, chdir to it; then, on context exit, add everything in it to git-annex"""
        temp_branch = os.path.join('tmp_wtree', worktree_group_id, str(uuid.uuid4()))
        temp_worktree_dir = os.path.join('tmp', temp_branch)
        self.execute_git(['worktree', 'add', '-b', temp_branch, temp_worktree_dir, self.get_first_commit()])
        dot_git = os.path.join(temp_worktree_dir, '.git')
        save_dot_git = util.file.slurp_file(dot_git)
        try:
            with contextlib2.ExitStack() as exit_stack:
                if chdir:
                    exit_stack.enter_context(util.file.pushd_popd(temp_worktree_dir))
                yield temp_branch, temp_worktree_dir
        finally:
            os.remove(dot_git)
            util.file.dump_file(dot_git, save_dot_git)
            if not keep and not util.file.keep_tmp():
                try:
                    self.execute_git(['worktree', 'remove', temp_worktree_dir])
                    self.execute_git(['branch', '-D', temp_branch])
                except Exception:
                    _log.warning('Could not remove temp worktree %s', temp_worktree_dir)

    def _call_f(self, f, arg_tuple):
        with self._in_tmp_worktree() as branch:
            try:
                f(arg_tuple)
            except Exception:
                return None
            
    # def in_tmp_worktrees(self, f, arg_tuples):
    #     """For each arg_tuple, call f(arg_tuple) on a temp branch in a temp worktree, checked out to an empty tree;
    #     then commit, and merge the results.
    #     """
    #     arg_tuples = tuple(arg_tuples)
    #     with contextlib2.ExitStack() as stk, \
    #          concurrent.futures.ThreadPoolExecutor(max_workers=min(len(arg_tuples), util.misc.available_cpu_count())) \
    #          as executor:
    #         worktree_group_id = '{:06d}'.format(random.randint(0,999999))
    #         temp_branches = [stk.enter_context(self._in_tmp_worktree(workree_group_id=worktree_group_id)) for arg_tuple in arg_tuples]
            
    #         result = filter(None, list(executor.map(functools.partial(GitAnnexTool._call_f, self=self, f=f), arg_tuples)))

# *** git-annex commands

    @_add_now_arg
    def add(self, fname):
        """Add a file to git-annex"""
        _log.debug('CALL TO ADD %s; batch status = %s', fname, self._batched_cmds)
        self.execute_batch(['add'], batch_args=(fname,))
        _log.debug('RETURNED FROM CALL TO ADD %s; batch status = %s', fname, self._batched_cmds)

    def add_cwd(self):
        """Add files in current working directory to git-annex"""
        self.execute(['add'])

    def initremote(self, name, remote_type, **kw):
        """Initialize a git-annex special remote"""
        remote_attrs = dict(kw)
        if 'encryption' not in remote_attrs: remote_attrs['encryption'] = 'none'
        self.execute(['initremote', name, 'type='+remote_type]+['='.join((k,v)) for k, v in remote_attrs.items()])

    def move(self, fname, to_remote_name):
        """Move a file from local repo to a remote"""
        self.execute(['move', fname, '--to', to_remote_name])

    @classmethod
    def _fromkey_preproc(cls, batched_calls):
        """Aux routine for fromkey().  Eliminate duplicate attempts to set a given file to point to the same key."""
        new_batched_calls = []
        fnames_seen = {}
        for batched_call in batched_calls:
            key, fname = batched_call.batch_args
            if fname in fnames_seen:
                util.misc.chk(fnames_seen[fname] == key)
                continue
            fnames_seen[fname] = key
            if os.path.lexists(fname):
                util.misc.chk(cls.lookupkey(fname) == key)
                continue
            new_batched_calls.append(batched_call)
                
        return new_batched_calls

    @_add_now_arg
    def fromkey(self, key, fname):
        """Manually set up a symlink to a given key as a given file.  If the target filename exists, and is anything
        other than an existing link to `key`, raise an error.
        """
        self.execute_batch(['fromkey', '--force'], batch_args=(key, os.path.abspath(fname)), preproc=self._fromkey_preproc)

    @classmethod
    def _get_link_into_annex(cls, f):
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

    @classmethod
    def is_link_into_annex(cls, f):
        """Test if `f` is a file controlled by git-annex (not necessarily present in local repo)."""
        return os.path.lexists(f) and not os.path.isdir(f) and cls._get_link_into_annex(f)[1]

    @classmethod
    def is_file_in_annex(cls, f):
        """Tests if `f` (possibly at the end of a chain of symlinks) is a file present in the local annex."""
        return cls.is_link_into_annex(f) and os.path.isfile(f)

    @classmethod
    def lookupkey(cls, f):
        """Get the git-annex key of an annexed file.  Note that, unlike git-annex-lookupkey, this looks at the file in the
        working copy, not in the index; and does not deal with unlocked files."""
        
        link_into_annex, target_of_link_into_annex = cls._get_link_into_annex(f)
        util.misc.chk(target_of_link_into_annex, '{} is not a link into annex'.format(os.path.abspath(f)))
        return os.path.basename(target_of_link_into_annex)

    def examinekey(self, key):
        """Return a dict of info that can be gleaned from the key itself.
        Most accurate would be to run 'git annex examinekey', but for speed we do the
        parsing ourselves.
        """
        key_attrs = collections.OrderedDict()
        key_attrs_str, key_name = key.rsplit('--', 1)
        key_attrs['key_name'] = key_name
        key_attrs_parts = key_attrs_str.split('-')
        key_attrs['backend'] = key_attrs_parts[0]
        key_attrs_size_part = [p for p in key_attrs_parts if p.startswith('s')]
        if key_attrs_size_part:
            key_attrs['size'] = int(key_attrs_size_part[0][1:])
        if key_attrs['backend'].startswith('MD5'):
            key_attrs['md5'] = key_name[:32]
        return key_attrs

    @staticmethod
    def _get_file_exts_for_key(fname, max_extension_length=5):
        """Determine the suffix of `fname` that would be included in the git-annex MD5E key, after the md5."""
        exts = ''
        num_exts = 0
        while True:
            fname, ext = os.path.splitext(fname)
            if not ext or len(ext) > (max_extension_length+1) or not ext[1:].isalnum():
                break
            exts = ext + exts
            num_exts += 1
            if num_exts >= 2:
                break
        return exts

    @staticmethod
    def construct_key(key_attrs, max_extension_length=5):
        """Construct a git-annex key from its attributes.  Currently only MD5E keys are supported.
        
        Args:
          key_attrs: map from key attribute name to value; must include mappings for the following keys -
             backend (currently must be MD5E), size, md5, fname.
        """
        key_attrs = copy.copy(key_attrs)
        key_attrs['md5'] = key_attrs['md5'].lower()
        util.misc.chk(key_attrs['backend'] == 'MD5E')
        _log.debug('constuct_key from %s', key_attrs)
        return '{backend}-s{size}--{md5}{exts}'.format(exts=GitAnnexTool._get_file_exts_for_key(key_attrs['fname'],
                                                                                                max_extension_length),
                                                       **key_attrs)

    def get_annexed_file_attrs(self, f):
        """Get annexed file info that can be determined quickly e.g. from the key"""
        raise NotImplemented()


    @_add_now_arg
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
                self.execute_batch(['get'], (os.path.basename(f),))

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

    @_add_now_arg
    def calckey(self, path, output_acceptor=None):
        """Compute the git-annex key for the given file"""
        our_output_acceptor = _ValHolder()
        self.execute(['calckey'], (path,), output_acceptor=output_acceptor or our_output_acceptor, now=not output_acceptor)
        return our_output_acceptor.val

    @_add_now_arg
    def registerurl(self, key, url):
        """Tell git-annex that `key` can be fetched from `url`"""
        self.execute_batch(['registerurl'], (key, url))

    @_add_now_arg
    def setpresentkey(self, key, remote_uuid, present):
        """Tell git-annex that `key` can be fetched from `url`"""
        self.execute_batch(['setpresentkey'], (key, remote_uuid, '1' if present else '0'))

    @_add_now_arg
    def register_key_in_remote_at_url(self, key, url, remote_uuid):
        """Record that `key` is present in remote `remote_uuid` at `url`"""
        self.setpresentkey(key, remote_uuid, True, now=False)
        self.registerurl(key, url, now=False)

    @_add_now_arg
    def checkpresentkey(self, key, output_acceptor=None, priority=-100):
        """Check if key is present in some remote"""
        our_output_acceptor = _ValHolder()
        self.execute(['checkpresentkey'], (key,), output_acceptor=output_acceptor or our_output_acceptor, now=not output_acceptor)
        return our_output_acceptor.val

    def _get_future(self, f, *args, **kw):
        if not self._executor:
            return f(*args, **kw)
        else:
            return self._executor.submit(f, *args, **kw)

    @staticmethod
    def _construct_key_for_file(file_path, key_attrs):
        key_attrs['md5'] = util.file.md5_for_file(file_path)
        return GitAnnexTool.construct_key(key_attrs=key_attrs)

    class Remote(object):
        
        """Abstract base cass for one instance of a particular git-annex special remote.
        """

        def __init__(self, remote_name, remote_uuid):
            self.remote_name = remote_name
            self.remote_uuid = remote_uuid

        def get_remote_uuid():
            """Return the git-annex uuid for this remote"""
            return self.remote_uuid

        def handles_url(self, url):
            """Return True if this remote handles this URL"""
            raise NotImplemented()

        def canonicalize_url(self, url, **kw):
            """Return a standard/canonical form of a given url.  For example, might convert a local
            pathname to file:///absolute/path.  The default implementation just returns the `url`.
            """
            return url

        def gather_filestats(self, ga_tool):
            """Gather filestats for URLs"""
            raise NotImplemented()

        @classmethod
        def load_external_special_remotes(cls, ga_tool, git_config):
            """Load external special remotes"""
            result = []
            for git_cfg_key, git_cfg_val in git_config.items():
                if git_cfg_key.startswith('remote.') and git_cfg_key.endswith('.annex-externaltype') \
                   and git_cfg_val==cls.EXTERNALTYPE:
                    remote_name = git_cfg_key.split('.')[1]
                    remote_uuid = git_config['remote.' + remote_name + '.annex-uuid']
                    result.append(cls(remote_name=remote_name, remote_uuid=remote_uuid))
            return result
        # end: def load_external_special_remotes(cls, ga_tool, git_config)

    # end: class Remote(object)

    class LocalDirRemote(Remote):
        """Files stored in local dir"""

        EXTERNALTYPE = 'ldir'

        def __init__(self, remote_name, remote_uuid):
            super(GitAnnexTool.LocalDirRemote, self).__init__(remote_name, remote_uuid)

        def handles_url(self, url):
            """Return True if this remote handles this URL"""
            return (uritools.isuri(url) and uritools.urisplit(url).scheme == 'file' and uritools.urisplit(url).authority == '') \
                or uritools.isabspath(url)

        def canonicalize_url(self, url, **kw):
            """Return a canonical form of a URL that this remote handles"""
            util.misc.chk(self.handles_url(url))
            if uritools.isuri(url):
                return url
            return uritools.uricompose(scheme='file', authority='', path=url)

        def gather_filestats(self, ga_tool):
            """Gather filestats for URLs that point to local files.

            For files that are git-annex links, get the git-annex key from the link, even if the file
            is not present in the local annex.

            For files for which we already have the md5, but not the size, quickly get the size and then compute
            the MD5E key.
            """

            _log.debug('LOCALDIRRMOTE gather')

            repo_root = ga_tool.get_repo_root()
            urls_in_this_repo = set()

            urls_with_new_keys = set()
            with ga_tool.batching() as ga_tool_calckey:
                for url, filestat in ga_tool._url2filestat.items():
                    _log.debug('LOOKING AT URL {} {}'.format(url, filestat))
                    if 'git_annex_key' in filestat: continue
                    if not self.handles_url(url): continue
                    canon_url = self.canonicalize_url(url)
                    
                    file_path = uritools.urisplit(canon_url).path
                    util.misc.chk(os.path.isabs(file_path))
                    
                    is_link_into_annex = ga_tool.is_link_into_annex(file_path)
                    if not (is_link_into_annex or os.path.isfile(file_path)): continue

                    urls_with_new_keys.add((url, canon_url))
                    
                    if is_link_into_annex:
                        filestat['git_annex_key'] = ga_tool.lookupkey(file_path)
                        if file_path.startswith(repo_root):
                            urls_in_this_repo.add(url)
                        continue

                    filestat['size'] = os.path.getsize(file_path)
                    if 'md5' not in filestat:
                        #filestat['md5'] = util.file.md5_for_file(file_path)
                        filestat['git_annex_key'] = ga_tool._get_future(ga_tool._construct_key_for_file, file_path,
                                                                        key_attrs=dict(backend='MD5E',
                                                                                       fname=os.path.basename(file_path),
                                                                                       **filestat))
                        continue
                    if 'md5' in filestat:
                        filestat['git_annex_key'] = ga_tool.construct_key(key_attrs=dict(backend='MD5E',
                                                                                         fname=os.path.basename(file_path),
                                                                                         **filestat))
                        continue

                    _log.debug('CALLING CALCKEY: %s', file_path)
                    ga_tool_calckey.calckey(file_path, output_acceptor=functools.partial(operator.setitem, filestat, 'git_annex_key'),
                                            now=False)
                # end: for url, filestat in ga_tool._url2filestat.items()
            # end: with ga_tool.batching() as ga_tool_calckey
            for url, canon_url in (urls_with_new_keys - urls_in_this_repo):
                ga_tool.register_key_in_remote_at_url(key=ga_tool._url2filestat[url]['git_annex_key'],
                                                      url=canon_url, remote_uuid=self.remote_uuid, now=False)
            
        # end: def gather_filestats(ga_tool):
    # end: class LocalDirRemote(object)

    class GsUriRemote(Remote):
        """Files stored in the google cloud"""

        EXTERNALTYPE = 'gs_uri'

        def __init__(self, remote_name, remote_uuid):
            super(GitAnnexTool.GsUriRemote, self).__init__(remote_name, remote_uuid)
            self.gcloud_tool = tools.gcloud.GCloudTool()

        def handles_url(self, url):
            """Return True if this remote handles this URL"""
            return uritools.isuri(url) and uritools.urisplit(url).scheme == 'gs'

        def gather_filestats(self, ga_tool):
            """Gather filestats for gs:// URLs.

            For files for which we already have the md5, but not the size, quickly get the size using gsutil ls.
            """

            gs_urls_needing_metadata = set()
            for url, filestat in ga_tool._url2filestat.items():
                if 'git_annex_key' in filestat: continue
                if not self.handles_url(url): continue
                gs_urls_needing_metadata.add(url)

            _log.debug('GS URLS NEEDING METADATA: %s', '\n'.join(map(str, gs_urls_needing_metadata)))
            url2attrs = self.gcloud_tool.get_metadata_for_objects(gs_urls_needing_metadata, ignore_errors=True)
            _log.debug('url2ATTRS=%s %s', url2attrs, ga_tool._url2filestat)

            for url in gs_urls_needing_metadata:
                filestat = ga_tool._url2filestat[url]
                if url in url2attrs:
                    filestat.update(size=url2attrs[url]['size'], md5=url2attrs[url]['md5'].lower())
                    filestat['git_annex_key'] = ga_tool.construct_key(key_attrs=dict(backend='MD5E',
                                                                                     size=filestat['size'],
                                                                                     md5=filestat['md5'],
                                                                                     fname=url))
                    ga_tool.register_key_in_remote_at_url(key=filestat['git_annex_key'],
                                                          url=url, remote_uuid=self.remote_uuid, now=False)

        # end: def gather_filestats(ga_tool):
    # end: class GsUriRemote(object)

    class DxRemote(Remote):
        """Files stored in DNAnexus"""

        EXTERNALTYPE = 'dnanexus'

        def __init__(self, remote_name, remote_uuid):
            super(GitAnnexTool.DxRemote, self).__init__(remote_name, remote_uuid)
            self.dx_tool = tools.dnanexus.DxTool()

        def handles_url(self, url):
            """Return True if this remote handles this URL"""
            return uritools.isuri(url) and uritools.urisplit(url).scheme == 'dx'

        @staticmethod
        def _key_from_dx_url(url):
            dxid = uritools.urisplit(url).authority
            dxf = dxpy.DXFile(dxid)
            descr = dxf.describe(fields={'parts'}, default_fields=True)
            if len(descr['parts']) == 1:
                md5 = descr['parts']['1']['md5']
            else:
                # TODO: add the use of a dxid-to-key cache in git (or more generally url-to-key?)
                with util.file.tempfname(suffix='_dx_md5') as t_file:
                    dxpy.dxfile_functions.download_dxfile(dxid, t_file)
                    md5 = util.file.md5_for_file(t_file)
            return GitAnnexTool.construct_key(key_attrs=dict(backend='MD5E',
                                                             size=descr['size'],
                                                             md5=md5,
                                                             fname=descr['name']))

        def gather_filestats(self, ga_tool):
            """Gather filestats for dx:// URLs.
            """

            dx_urls_needing_metadata = set()
            for url, filestat in ga_tool._url2filestat.items():
                if 'git_annex_key' in filestat: continue
                if not self.handles_url(url): continue
                dx_urls_needing_metadata.add(url)

            _log.debug('DX URLS NEEDING METADATA: %s', '\n'.join(map(str, dx_urls_needing_metadata)))

            for url in dx_urls_needing_metadata:
                filestat = ga_tool._url2filestat[url]
                filestat['git_annex_key'] = ga_tool._get_future(self._key_from_dx_url, url)
                ga_tool.register_key_in_remote_at_url(key=filestat['git_annex_key'],
                                                      url=url, remote_uuid=self.remote_uuid, now=False)

        # end: def gather_filestats(ga_tool):
    # end: class DxRemote(object)

    REMOTE_TYPES = [LocalDirRemote, GsUriRemote, DxRemote]

    def get_remotes(self):
        """Load remotes from the git configuration"""
        git_config = self.get_git_config()
        return functools.reduce(operator.concat, [list(remote_type.load_external_special_remotes(self, git_config))
                                                  for remote_type in self.REMOTE_TYPES], [])

    def get_key_for_url(self, url, error_if_missing=True):
        """Return git-annex key for the given URL.  If not known, raise an exception if error_if_missing is True,
        else return None.
        """
        key = self._url2filestat.get(url, {}).get('git_annex_key', None)
        if not key and error_if_missing:
            raise RuntimeError('No key for url {}'.format(url))
        return key

    def get_url2filestat(self):
        """Return the internal mapping of URLs to filestats"""
        return self._url2filestat
    
    @staticmethod
    def _chk_present(val, key):
        key = util.misc.maybe_wait_for_result(key)
        util.misc.chk_eq(val,  '1', message='failed to import key {}'.format(key))

    @_add_now_arg
    def import_urls(self, urls, ignore_non_urls=False, check_imported=True):
        """Imports `urls` into git-annex.

        Args:
          urls: iterable of urls, which are strings denoting either local or cloud files.
          ignore_non_urls: if True, urls we don't recognize as urls will be silently ignored; if False, they'll trigger an error.
          check_imported: double-check that each URL has been successfully imported
          now: if False, effects may be delayed until batching context ends

        Effects:
          each valid url in `urls` becomes known to git-annex, i.e. the presence of its key in the repository designated by
          the url is recorded in the git-annex branch.

        Returns:
          map from url to filestat.  For each valid url in `urls`, the resulting filestat will contain at least the following
          keys: 'git_annex_key' .
        """

        urls = sorted(set(urls))

        for url in urls:
            if url not in self._url2filestat:
                self._url2filestat[url] = collections.OrderedDict()

        remotes = self.get_remotes()

        for remote in remotes:
            remote.gather_filestats(ga_tool=self)
        # end: for remote in remotes

        _log.debug('GOT RESULTS: %s', {k: v for k, v in self._url2filestat.items() if v})

        for url in urls:
            util.misc.chk(ignore_non_urls  or  ('git_annex_key' in self._url2filestat[url]),
                          'no key for {}: {}'.format(url, self._url2filestat))
            if False and check_imported and 'git_annex_key' in self._url2filestat[url]:
                key = self._url2filestat[url]['git_annex_key']
                self.checkpresentkey(key=key,
                                     output_acceptor=functools.partial(GitAnnexTool._chk_present, key=key),
                                     now=False)
    # end: def import_urls(self, urls, ignore_non_urls=False, check_imported=True):

    @_add_now_arg
    def from_url(self, url, fname):
        """Import url and set fname to point to the resulting key"""
        self.import_urls([url])
        self.fromkey(self._url2filestat[url], os.path.abspath(fname))

# end: class GitAnnexTool
