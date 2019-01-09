'''
    Tool wrapper for git-annex
'''

# * imports

import logging
import collections
import os
import os.path
import stat
import subprocess
import shutil
import random
import shlex
import tempfile
import pipes
import time
import functools
import copy
import warnings

import contextlib2 as contextlib

import tools
import util.file
import util.misc
import util.version

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
class _ValHolder(object):
    def __init__(self, val=None):
        self.val = val

_BatchedCall = collections.namedtuple('_BatchedCall', ['batch_args', 'output_holder'])

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
        self._batched_cmds = None
        _log.info('Created GitAnnexTool %s', id(self))

    def version(self):
        return TOOL_VERSION

    def _get_bin_dir(self):
        """Return the directory containing the git-annex binary"""
        bin_dir = os.path.dirname(self.install_and_get_path())
        util.misc.chk(self.is_installed())
        return bin_dir

    def _get_run_env(self):
        """Return the environment in which git-annex should be run.  This environment will include in its PATH
        the binaries for git-annex external special remote implementations that we add."""
        if not hasattr(self, '_run_env'):
            # ensure that the git and git-annex binaries we installed are first in PATH.
            self._run_env = dict(os.environ, PATH=':'.join((self._get_bin_dir(),
                                                            os.path.join(util.version.get_project_path(),
                                                                         'tools', 'git-annex-remotes'),
                                                            os.environ['PATH'])))
        return self._run_env

    def execute(self, *args, **kwargs):
        """Run the given git or git-annex command.  If batching is in effect (see self.batching()), the command is stored and will be
        executed when the batching context exits.
        
        Args:
            args: list of arguments too the command
            tool: git-annex or git, defaults to git-annex
            batch_args: if given, tuple of args to a batch git-annex command
            output_holder: if given, the output of the command will be stored into this _ValHolder object.
              Note that if batching is in effect, this will only happen after the batching context closes.
        """
        with contextlib.ExitStack() as stack:
            if self._batched_cmds is None:
                self = stack.enter_context(self.batching())
            self._add_command_to_batch(*args, **kwargs)

    def _add_command_to_batch(self, args, tool='git-annex', batch_args=None, output_holder=None):
        """Add command to current batch"""
        tool_cmd = (os.path.join(self._get_bin_dir(), tool),) + tuple(map(str, args))
        self._batched_cmds[tool_cmd].append(_BatchedCall(batch_args=batch_args, output_holder=output_holder))

    def _execute_batched_commands(self):
        """Run any saved batched commands.

        Note that if one command fails, others won't be run.
        """

        # TODO: figure out correct behavior if some cmds fail

        _log.debug('RUNNING BATCH CMDS: %s', self._batched_cmds)
        while self._batched_cmds:
            tool_cmd, batch_calls = self._batched_cmds.popitem()
            with contextlib.ExitStack() as stack:
                subprocess_call_args = {}
                if batch_calls[0].batch_args:
                    _log.info('calling batch calls: %s %s', tool_cmd, batch_calls)
                    batch_inps_fname = stack.enter_context(util.file.tempfname(prefix='git-annex-batch-inputs'))
                    util.file.dump_file(batch_inps_fname,
                                        '\n'.join([' '.join(map(str, batch_call.batch_args)) for batch_call in batch_calls]))
                    batch_inps_file = stack.enter_context(open(batch_inps_fname))
                    subprocess_call_args.update(stdin=batch_inps_file)

                subprocess_func = getattr(subprocess, 'check_output' if batch_calls[0].output_holder else 'check_call')
                result = subprocess_func(tool_cmd, env=self._get_run_env(), **subprocess_call_args)

                if batch_calls[0].output_holder:
                    if len(batch_calls) == 1:
                        batch_calls[0].output_holder.val = result
                    else:
                        result_lines = result.strip().split()
                        util.misc.chk(len(result_lines) == len(batch_calls))
                        for batch_call, result_line in zip(batch_calls, result_lines):
                            batch_call.output_holder.val = result_line
                
    @contextlib.contextmanager
    def batching(self):
        """Create a batching context.  Commands issued within the batching context _may_ be delayed and issued in batches
        for efficiency; we guarantee that when the batching context exits, the commands have been executed."""
        self = copy.copy(self)
        self._batched_cmds = collections.defaultdict(list)
        yield self
        if not self._batched_cmds:
            warnings.warn("Empty batching context -- did you forget to use the yielded GitAnnexTool()?", RuntimeWarning)
        self._execute_batched_commands()

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

    def initremote_external(self, remote_name, externaltype):
        """Init an external special remote ldir"""
        self.execute(['initremote', remote_name, 'type=external', 'externaltype={}'.format(externaltype), 'encryption=none'])

    def add(self, fname):
        """Add a file to git-annex"""
        _log.debug('CALL TO ADD %s; batch status = %s', fname, self._batched_cmds)
        self.execute(['add', '--batch'], batch_args=(fname,))
        _log.debug('RETURNED FROM CALL TO ADD %s; batch status = %s', fname, self._batched_cmds)

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
        while True:
            fname, ext = os.path.splitext(fname)
            if not ext or len(ext) > (max_extension_length+1):
                break
            exts = ext + exts
        return exts

    @staticmethod
    def construct_key(key_attrs, max_extension_length=5):
        """Construct a git-annex key from its attributes.  Currently only MD5E keys are supported."""
        util.misc.chk(key_attrs['backend'] == 'MD5E')
        _log.info('constuct_key from %s', key_attrs)
        return '{backend}-s{size}--{md5}{exts}'.format(exts=GitAnnexTool._get_file_exts_for_key(key_attrs['fname'],
                                                                                                max_extension_length),
                                                       **key_attrs)

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

    def calckey(self, path, output_holder):
        """Compute the git-annex key for the given file"""
        self.execute(['calckey', '--batch'])
        

    def import_urls(self, urls):
        """Imports `urls` into git-annex.  urls can be local or cloud paths.

        Returns:
           map from url to git-annex key for the contents of the URL.
        """
        
        pass

    def _gather_fmdata_from_local_files(lcpath2fmdata):
        """Gather file metadata from local files."""
        for file_path, fmdata in lcpath2fmdata.items():
            if not file_path.startswith('/') or all(f in fmdata for f in ('md5', 'size')) or \
               not os.path.isfile(file_path):
                continue
            fmdata['size'] = os.path.getsize(file_path)
            fmdata['md5'] = _md5_for_file(file_path)

        

# end: class GitAnnexTool
