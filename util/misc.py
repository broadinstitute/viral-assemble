'''A few miscellaneous tools. '''
from __future__ import print_function, division  # Division of integers with / should never round!
import collections
import contextlib
import itertools, functools, operator
import logging
import os, os.path
import re
import subprocess
import threading
import multiprocessing
import sys
import copy
import yaml, json
import inspect
import time
import builtins
import argparse
import concurrent.futures
import json
import signal
import socket

import util.file

import sqlitedict
import psutil

log = logging.getLogger(__name__)

__author__ = "dpark@broadinstitute.org"

MAX_INT32 = (2 ** 31)-1

@contextlib.contextmanager
def timer(prefix):
    start = time.time()
    yield
    finish = time.time()
    elapsed = '{:.2f}'.format(finish - start)
    print(prefix + ' - ' + elapsed, file=sys.stderr)

def get_module_name(f):
    """Get module name of a function"""
    return os.path.splitext(os.path.relpath(inspect.getsourcefile(util.misc.unwrap(f)),
                                            util.version.get_project_path()))[0].replace(os.sep, '.')

def memoize(obj):
    """Function decorator that caches function return values for given arguments, and
    returns results from the cache when available. Must not be used for functions that are non-deterministic, have
    side effects, or depend on external state that may change.

    See also: memoize_persist()
    """

    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = "".join([str(args),str(kwargs)])
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer

def identity(x):
    """The identity function"""
    return x

_caches = []
def _commit_caches():
    for cache in _caches:
        cache.commit()

def cleanup_misc():
    """Do cleanup tasks for the `util.misc` module"""
    _commit_caches()

def memoize_persist(to_picklable=identity, from_picklable=identity):
    """Function decorator that caches function return values for given arguments, and
    returns results from the cache when available.  The cache can persist on disk across
    command executions.  Must not be used for functions that are non-deterministic, have
    side effects, or depend on external state that may change.

    Note: if the implementation of the cached function or any function it calls changes,
    the cache must be manually cleared!
    """

    def memoize_persist_do(orig_fn):
        cache_dir = os.path.join(util.file.get_cache_dir(), 'vngs', get_module_name(orig_fn))
        util.file.mkdir_p(cache_dir)

        # We store arg-to-result mappings using sqlitedict .
        # To reduce cache contention, we split the database into many sqlite files based on key.

        def _get_cache(key, key_hash2cache = {}):
            """Return the SqliteDict for the cache responsible for storing result of `key`."""
            key_hash = '{:03d}'.format(hash(key) % 1000)

            if key_hash not in key_hash2cache:
                db_dir = os.path.join(cache_dir, orig_fn.__name__, key_hash[0], key_hash[1])
                util.file.mkdir_p(db_dir)
                cache = sqlitedict.SqliteDict(os.path.join(db_dir, key_hash[2]))
                _caches.append(cache)
                key_hash2cache[key_hash] = cache
            return key_hash2cache.get[key_hash]

        @wraps(orig_fn)
        def memoizer(*args, **kwargs):
            key = "".join([str(args),str(sorted(kwargs.items()))])
            try:
                cache = _get_cache(key)
            except Exception:
                cache = {}
            if key in cache:
                return from_picklable(cache[key])
            result = orig_fn(*args, **kwargs)
            try:
                cache[key] = to_picklable(result)
            except Exception:
                log.warning('memoize_persist: could not cache result of %s for args=%s kwargs=%s',
                            orig_fn.__name__, args, kwargs)
                pass
            return result
        return memoizer
    return memoize_persist_do
# end: def memoize_persist(to_picklable=identity, from_picklable=identity)

def unique(items):
    ''' Return unique items in the same order as seen in the input. '''
    seen = set()
    for i in items:
        if i not in seen:
            seen.add(i)
            yield i


def histogram(items):
    ''' I count the number of times I see stuff and return a dict of counts. '''
    out = {}
    for i in items:
        out.setdefault(i, 0)
        out[i] += 1
    return out


def freqs(items, zero_checks=None):
    ''' Given a list of comparable, non-unique items, produce an iterator of
            (item, count, freq) tuples.
            item is a unique instance of one of the items seen on input
            count is a positive integer describing the number of times this item was observed
            freq is count / sum(count) for all outputs.
        If zero_checks is specified, then the output iterator will emit tuples for the
        items in zero_checks even if they do not appear in the input.  If they are not in
        the input, they will be emitted with a zero count and freq.
        See histogram(items)
    '''
    zero_checks = zero_checks or set()

    tot = 0
    out = {}
    for i in items:
        out.setdefault(i, 0)
        out[i] += 1
        tot += 1
    for k, v in out.items():
        yield (k, v, float(v) / tot)
    for i in zero_checks:
        if i not in out:
            yield (i, 0, 0.0)


def intervals(i, n, l):
    ''' Divide something of length l into n equally sized parts and return the
        start-stop values of the i'th part.  Values are 1-based.  Each part
        will be adjacent and non-overlapping with the next part. i must be a
        number from 1 to n.
    '''
    assert 1 <= i <= n and l >= n
    part_size = l // n
    start = 1 + part_size * (i - 1)
    stop = part_size * i
    if i == n:
        stop = l
    return (start, stop)

# from http://stackoverflow.com/a/312467


def pairwise(iterable):
    """ from itertools recipes
        s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = itertools.tee(iterable)
    next(b, None)
    if hasattr(itertools, 'izip'):
        # Python 2
        return itertools.izip(a, b)
    else:
        # Python 3
        return zip(a, b)


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    it = iter(iterator)
    item = list(itertools.islice(it, batch_size))
    while item:
        yield item
        item = list(itertools.islice(it, batch_size))


def list_contains(sublist, list_):
    """Tests whether sublist is contained in full list_."""

    for i in range(len(list_)-len(sublist)+1):
        if sublist == list_[i:i+len(sublist)]:
            return True
    return False


try:
    from subprocess import run
except ImportError:
    CompletedProcess = collections.namedtuple(
        'CompletedProcess', ['args', 'returncode', 'stdout', 'stderr'])

    def run(args, stdin=None, stdout=None, stderr=None, shell=False,
            env=None, cwd=None, timeout=None, check=False, executable=None):
        '''A poor man's substitute of python 3.5's subprocess.run().

        Should only be used for capturing stdout. If stdout is unneeded, a
        simple subprocess.call should suffice.
        '''
        assert stdout is not None, (
            'Why are you using this util function if not capturing stdout?')

        stdout_pipe = stdout == subprocess.PIPE
        stderr_pipe = stderr == subprocess.PIPE
        # A little optimization when we don't need temporary files.
        if stdout_pipe and (
                stderr == subprocess.STDOUT or stderr is None):
            try:
                output = subprocess.check_output(
                    args, stdin=stdin, stderr=stderr, shell=shell,
                    env=env, cwd=cwd, executable=executable)
                return CompletedProcess(args, 0, output, b'')
            # Py3.4 doesn't have stderr attribute
            except subprocess.CalledProcessError as e:
                if check:
                    raise
                returncode = e.returncode
                stderr_text = getattr(e, 'stderr', b'')
                return CompletedProcess(args, e.returncode, e.output, stderr_text)

        # Otherwise use temporary files as buffers, since subprocess.call
        # cannot use PIPE.
        if stdout_pipe:
            stdout_fn = util.file.mkstempfname('.stdout')
            stdout = open(stdout_fn, 'wb')
        if stderr_pipe:
            stderr_fn = util.file.mkstempfname('.stderr')
            stderr = open(stderr_fn, 'wb')
        try:
            returncode = subprocess.call(
                args, stdin=stdin, stdout=stdout,
                stderr=stderr, shell=shell, env=env, cwd=cwd,
                executable=executable)
            if stdout_pipe:
                stdout.close()
                with open(stdout_fn, 'rb') as f:
                    output = f.read()
            else:
                output = ''
            if stderr_pipe:
                stderr.close()
                with open(stderr_fn, 'rb') as f:
                    error = f.read()
            else:
                error = ''
            if check and returncode != 0:
                print(output.decode("utf-8"))
                print(error.decode("utf-8"))
                try:
                    raise subprocess.CalledProcessError(
                        returncode, args, output, error) #pylint: disable-msg=E1121
                except TypeError: # py2 CalledProcessError does not accept error
                    raise subprocess.CalledProcessError(
                        returncode, args, output)
            return CompletedProcess(args, returncode, output, error)
        finally:
            if stdout_pipe:
                stdout.close()
                os.remove(stdout_fn)
            if stderr_pipe:
                stderr.close()
                os.remove(stderr_fn)


def run_and_print(args, stdout=None, stderr=subprocess.STDOUT,
                  stdin=None, shell=False, env=None, cwd=None,
                  timeout=None, silent=False, buffered=False, check=False,
                  loglevel=None):
    '''Capture stdout+stderr and print.

    This is useful for nose, which has difficulty capturing stdout of
    subprocess invocations.
    '''
    if loglevel:
        silent = False
    if not buffered:
        if check and not silent:
            try:
                result = run(
                    args,
                    stdin=stdin,
                    stdout=subprocess.PIPE,
                    stderr=stderr,
                    env=env,
                    cwd=cwd,
                    timeout=timeout,
                    check=check
                )
                print(result.stdout.decode('utf-8'))
                try:
                    print(result.stderr.decode('utf-8'))
                except AttributeError:
                    pass
            except subprocess.CalledProcessError as e:
                if loglevel:
                    try:
                        log.log(loglevel, result.stdout.decode('utf-8'))
                    except NameError:
                        # in some situations, result does not get assigned anything
                        pass
                    except AttributeError:
                        log.log(loglevel, result.output.decode('utf-8'))
                else:
                    print(e.output.decode('utf-8'))
                    try:
                        print(e.stderr.decode('utf-8'))
                    except AttributeError:
                        pass
                    sys.stdout.flush()
                raise(e)
        else:
            result = run(args, stdin=stdin, stdout=subprocess.PIPE,
                         stderr=stderr, env=env, cwd=cwd,
                         timeout=timeout, check=check)
            if not silent and not loglevel:
                print(result.stdout.decode('utf-8'))
                sys.stdout.flush()
            elif loglevel:
                log.log(loglevel, result.stdout.decode('utf-8'))

    else:
        CompletedProcess = collections.namedtuple(
        'CompletedProcess', ['args', 'returncode', 'stdout', 'stderr'])

        process = subprocess.Popen(args, stdin=stdin, stdout=subprocess.PIPE,
                                    stderr=stderr, env=env,
                                    cwd=cwd)
        output = []
        while process.poll() is None:
            for c in iter(process.stdout.read, b""):
                output.append(c)
                if not silent:
                   print(c.decode('utf-8'), end="") # print for py3 / p2 from __future__
                sys.stdout.flush() # flush buffer to stdout

        # in case there are still chars in the pipe buffer
        for c in iter(lambda: process.stdout.read(1), b""):
           if not silent:
               print(c, end="")
           sys.stdout.flush() # flush buffer to stdout

        if check and process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, args,
                                                b''.join(output))
        result = CompletedProcess(
            args, process.returncode, b''.join(output), b'')

    return result


def run_and_save(args, stdout=None, stdin=None,
                 outf=None, stderr=None, preexec_fn=None,
                 close_fds=False, shell=False, cwd=None, env=None, check=True):
    assert outf is not None

    sp = subprocess.Popen(args,
                          stdin=stdin,
                          stdout=outf,
                          stderr=subprocess.PIPE,
                          preexec_fn=preexec_fn,
                          close_fds=close_fds,
                          shell=shell,
                          cwd=cwd,
                          env=env)
    out, err = sp.communicate()

    # redirect stderror to stdout so it can be captured by nose
    if err:
        sys.stdout.write(err.decode("UTF-8"))

    if sp.returncode != 0 and check:
        raise subprocess.CalledProcessError(sp.returncode, str(args[0]))

    return sp


class FeatureSorter(object):
    ''' This class helps sort genomic features. It's not terribly optimized
        for speed or anything. Slightly inspired by calhoun's MultiSequenceRangeMap.
    '''
    def __init__(self, collection=None):
        self.seqids = []
        self.seq_to_features = {}
        self.seq_to_breakpoints = {}
        self.dirty = False
        if collection is not None:
            for args in collection:
                self.add(*args)

    def add(self, c, start, stop, strand='+', other=None):
        ''' Add a "feature", which is a chrom,start,stop,strand tuple (with
            optional other info attached)
        '''
        if c not in self.seq_to_features:
            self.seqids.append(c)
            self.seq_to_features[c] = []
            self.seq_to_breakpoints[c] = set()
            #self.seq_to_breakpoints[c].add(1) # do we want the empty interval in front?
        self.seq_to_features[c].append((int(start), int(stop), strand, other))
        self.seq_to_breakpoints[c].add(start)
        self.seq_to_breakpoints[c].add(stop+1)
        self.dirty = True

    def _cleanup(self):
        if self.dirty:
            self.dirty = False
            for c in self.seqids:
                self.seq_to_features[c].sort()

    def get_seqids(self):
        return tuple(self.seqids)

    def get_features(self, c=None, left=0, right=float('inf')):
        ''' Get all features on all chromosomes in sorted order. Chromosomes
            are emitted in order of first appearance (via add). Features on
            each chromosome are emitted in start, then stop order.  If
            boundaries are specified, we restrict to features that contain
            the specified interval.
        '''
        self._cleanup()
        if c is not None:
            seqlist = [c]
        else:
            seqlist = self.seqids
        for c in seqlist:
            for start, stop, strand, other in self.seq_to_features[c]:
                if stop>=left and start<=right:
                    yield (c, start, stop, strand, other)

    def get_intervals(self, c=None):
        ''' Get all intervals on the reference where the overlapping feature
            set remains the same. Output will be sorted, adjacent intervals
            and will describe how many and which features overlap it.
        '''
        self._cleanup()
        if c is not None:
            seqlist = [c]
        else:
            seqlist = self.seqids
        for c in seqlist:
            for left, right in pairwise(sorted(self.seq_to_breakpoints[c])):
                right = right - 1
                features = list(self.get_features(c, left, right))
                yield (c, left, right, len(features), features)


def available_cpu_count():
    """
    Return the number of available virtual or physical CPUs on this system.
    The number of available CPUs can be smaller than the total number of CPUs
    when the cpuset(7) mechanism is in use, as is the case on some cluster
    systems.

    Adapted from http://stackoverflow.com/a/1006301/715090
    """

    cgroup_cpus = MAX_INT32
    try:
        def get_cpu_val(name):
            return float(util.file.slurp_file('/sys/fs/cgroup/cpu/cpu.'+name).strip())
        cfs_quota = get_cpu_val('cfs_quota_us')
        if cfs_quota > 0:
            cfs_period = get_cpu_val('cfs_period_us')
            log.debug('cfs_quota %s, cfs_period %s', cfs_quota, cfs_period)
            cgroup_cpus = max(1, int(cfs_quota / cfs_period))
    except Exception as e:
        pass

    proc_cpus = MAX_INT32
    try:
        with open('/proc/self/status') as f:
            status = f.read()
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$', status)
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                proc_cpus = res
    except IOError:
        pass

    log.debug('cgroup_cpus %d, proc_cpus %d, multiprocessing cpus %d',
              cgroup_cpus, proc_cpus, multiprocessing.cpu_count())
    return min(cgroup_cpus, proc_cpus, multiprocessing.cpu_count())

def sanitize_thread_count(threads=None, tool_max_cores_value=available_cpu_count):
    ''' Given a user specified thread count, this function will:
        - ensure that 1 <= threads <= available_cpu_count()
        - interpret None values to mean max available cpus
            unless PYTEST_XDIST_WORKER_COUNT is defined as an environment
            variable, in which case we always return 1
        - allow for various, tool-specific ways of specifying
            max available cores (tool_max_value)
        tool_max_cores_value can be one of:
            available_cpu_count - this function will return available_cpu_count()
            any other value - this function will return that value.
            some commonly used values for tools are -1, 0, and None.
    '''
    if 'PYTEST_XDIST_WORKER_COUNT' in os.environ:
        threads = 1

    max_cores = available_cpu_count()

    if threads is None:
        threads = max_cores

    assert type(threads) == int

    if threads >= max_cores:
        if tool_max_cores_value == available_cpu_count:
            threads = max_cores
        else:
            threads = tool_max_cores_value
    else:
        if threads < 1:
            threads = 1
    return threads

def which(application_binary_name):
    """
        Similar to the *nix "which" command,
        this function finds the first executable binary present
        in the system PATH for the binary specified.
        It differs in that it resolves symlinks.
    """
    path=os.getenv('PATH')
    for path in path.split(os.path.pathsep):
        full_path=os.path.join(path, application_binary_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            link_resolved_path = os.path.realpath(full_path)
            return link_resolved_path

def is_nonatom_iterable(x, atom_types=(str,)):
    '''Tests whether `x` is an Iterable other than a string'''
    return isinstance(x, collections.Iterable) and not isinstance(x,atom_types)

_str_type = getattr(builtins, 'basestring', getattr(builtins, 'str'))

def make_seq(x, atom_types=(_str_type,)):
    '''Return a tuple containing the items in `x`, or containing just `x` if `x` is a non-atomic iterable.  Convenient
    for uniformly writing iterations over parameters that may be passed in as either an item or a tuple/list of items.
    Note that if `x` is an iterator, it will be concretized.  `atom_types` gives the type(s) to treat as atomic.'
    '''
    return tuple(x) if is_nonatom_iterable(x, atom_types) else (x,)

def load_yaml_or_json(fname):
    '''Load a dictionary from either a yaml or a json file'''
    with open(fname) as f:
        if fname.upper().endswith('.YAML'): return yaml.safe_load(f) or {}
        if fname.upper().endswith('.JSON'): return json.load(f) or {}
        raise TypeError('Unsupported dict file format: ' + fname)


def load_config(cfg, include_directive='include', std_includes=(), param_renamings=None):
    '''Load a configuration, with support for some extra functionality that lets project configurations evolve
    without breaking backwards compatibility.

    The configuration `cfg` is either a dict (possibly including nested dicts) or a yaml/json file containing one.
    A configuration parameter or config param is a sequence of one or more keys; the value of the corresponding
    parameter is accessed as "cfg[k1][k2]...[kN]".  Note, by "parameter" we denote the entire sequence of keys.

    This function implements extensions to the standard way of specifying configuration parameters via (possibly nested)
    dictionaries.  These extensions make it easier to add or rename config params without breaking backwards
    compatibility.

    One extension lets config files include other config files, and lets you specify "standard" config file(s) to
    be included before all others.  If the "default" config file from the project distribution is made a standard
    include, new parameters can be added to it (and referenced from project code) without breaking compatibility
    with old config files that omit these parameters.

    Another extension lets you, when loading a config file, recognize parameters specified under old or legacy names.
    This lets you change parameter names in new program versions while still accepting legacy config files that
    use older names.

    Args:
       cfg: either a config mapping, or the name of a file containing one (in yaml or json format).
         A config mapping is just a dict, possibly including nested dicts, specifying config params.
         The mapping can include an entry pointing to other config files to include.
         The key of the entry is `include_directive`, and the value is a filename or list of filenames of config files.
         Relative filenames are interpreted relative to the directory containing `cfg`, if `cfg` is a filename,
         else relative to the current directory.  Any files from `std_includes` are prepended to the list of
         included config files.  Parameter values from `cfg` override ones from any included files, and parameter values
         from included files listed later in the include list override parameter values from included files listed
         earlier.

       include_directive: key used to specify included config files
       std_includes: config file(s) implicitly included before all others and before `cfg`
       param_renamings: optional map of old/legacy config param names to new ones.  'Param names' here are
         either keys or sequences of keys.  Example value: {'trinity_kmer_size': ('de_novo_assembly', 'kmer_size')};
         new code can access the parameter as cfg["de_novo_assembly"]["kmer_size"] while legacy users can keep
         specifying it as "trinity_kmer_size: 31".
    '''

    param_renamings = param_renamings or {}

    result = dict()

    base_dir_for_includes = None
    if isinstance(cfg, str):
        cfg_fname = os.path.realpath(cfg)
        base_dir_for_includes = os.path.dirname(cfg_fname)
        cfg = load_yaml_or_json(cfg_fname)

    def _update_config(config, overwrite_config):
        """Recursively update dictionary config with overwrite_config.

        Adapted from snakemake.utils.
        See
        http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
        for details.

        Args:
          config (dict): dictionary to update
          overwrite_config (dict): dictionary whose items will overwrite those in config

        """

        def _update(d, u):
            def fix_None(v): return {} if v is None else v
            for (key, value) in u.items():
                if (isinstance(value, collections.Mapping)):
                    d[key] = _update(fix_None(d.get(key, {})), value)
                else:
                    d[key] = fix_None(value)
            return d

        _update(config, overwrite_config)


    includes = make_seq(std_includes) + make_seq(cfg.get(include_directive, []))
    for included_cfg_fname in includes:
        if (not os.path.isabs(included_cfg_fname)) and base_dir_for_includes:
            included_cfg_fname = os.path.join(base_dir_for_includes, included_cfg_fname)
        _update_config(result, load_config(cfg=included_cfg_fname, 
                                           include_directive=include_directive,
                                           param_renamings=param_renamings))

    # mappings in the current (top-level) config override any mappings from included configs
    _update_config(result, cfg)

    # load any params specified under legacy names, for backwards compatibility
    param_renamings_seq = dict(map(lambda kv: map(make_seq, kv), param_renamings.items()))

    for old_param, new_param in param_renamings_seq.items():

        # handle chains of param renamings
        while new_param in param_renamings_seq:
            new_param = param_renamings_seq[new_param]

        old_val = functools.reduce(lambda d, k: d.get(k, {}), old_param, result)
        new_val = functools.reduce(lambda d, k: d.get(k, {}), new_param, result)

        if old_val != {} and new_val == {}:
            _update_config(result, functools.reduce(lambda d, k: {k: d}, new_param[::-1], old_val))
            log.warning('Config param {} has been renamed to {}; old name accepted for now'.format(old_param, new_param))

    return result

def flatten(d, types_to_flatten=(tuple, list)):
    """Flatten a data structure `d` into a flat list.  Any elements of `d` that are instances of one of the iterable types in 
    `types_to_flatten` will be recursively flattened."""
    return [d] if not isinstance(d, types_to_flatten) else \
        functools.reduce(operator.concat, list(map(functools.partial(flatten, types_to_flatten=types_to_flatten), d)), [])

FullArgSpec = collections.namedtuple('FullArgSpec',
                                     'args, varargs, varkw, defaults, kwonlyargs, kwonlydefaults, annotations')
def getfullargspec(func):
    """Backport of inspect.getfullargspec() from Python 3 to Python 2"""
    if hasattr(inspect, 'getfullargspec'): return inspect.getfullargspec(func)
    argspec = inspect.getargspec(func)
    return FullArgSpec(argspec[0], argspec[1], argspec[2], argspec[3], kwonlyargs=[], kwonlydefaults={}, annotations={})

def getnamedargs(func):
    """Return a list of named args -- both positional and keyword-only -- of the given function."""
    argspec = getfullargspec(func)
    return flatten(argspec.args + argspec.kwonlyargs)

def as_type(val, types):
    """Try converting `val`to each of `types` in turn, returning the first one that succeeds."""
    errs = []
    for type in make_seq(types):
        try:
            return type(val)
        except Exception as e:
            errs.append(e)
            pass
    raise TypeError('Could not convert {} to any of {}: {}'.format(val, types, errs))

def subdict(d, keys):
    """Return a newly allocated shallow copy of a mapping `d` restricted to keys in `keys`."""
    d = collections.OrderedDict(d)
    keys = set(keys)
    return collections.OrderedDict([(k, v) for k, v in d.items() if k in keys])

def chk(condition, message='Check failed', exc=RuntimeError):
    """Check a condition, raise an exception if bool(condition)==False, else return `condition`."""
    if not condition:
        raise exc(message)
    return condition

def chk_eq(a, b, message='Check failed: {} != {}', exc=RuntimeError):
    """Check that a == b, raise an exception if not, else return True."""
    return chk(a == b, message.format(a, b), exc=exc)

def from_module(module, names):
    """Shorthand for assigning several functions from a module to local vars,
    to avoid code littered with full module paths."""
    return [getattr(module, name) for name in names.split()]

def first_non_None(*args):
    """Return the first of args that is not None; raise ValueError if all are None."""
    for arg in args:
        if arg is not None:
            return arg
    raise ValueError('No non-None args')

def maybe_decode(s):
    """Return a unicode version of `s`, decoding if needed.  For Python 2/3 compatibility."""
    return s.decode() if hasattr(s, 'decode') else s

def transform_json_data(val, node_handler, path=()):
    """Transform a parsed json structure, by replacing nodes for which `node_handler` returns
    an object different from `val`, with that object.  If `node_handler` has a named arg `path`,
    the path from the root will be passed via that arg (as a tuple).  If `node_handler` has a named
    arg `is_leaf`, it will set to True when node_handler is called for a leaf, to False otherwise.
    """
    recurse = functools.partial(transform_json_data, node_handler=node_handler)
    node_handler_named_args = getnamedargs(node_handler)
    node_handler_args = {}
    if 'path' in node_handler_named_args:
        node_handler_args['path'] = path
    is_leaf = not isinstance(val, (list, collections.Mapping))
    if 'is_leaf' in node_handler_named_args:
        node_handler_args['is_leaf'] = is_leaf
    handled_val = node_handler(val, **node_handler_args)
    if handled_val is not val:
        #log.debug('resolved %s to %s', val, handled_val)
        return handled_val
    elif isinstance(val, list):
        return [recurse(val=v, path=path+(i,)) for i, v in enumerate(val)]
    elif isinstance(val, collections.Mapping):
        return collections.OrderedDict([(k, recurse(val=v, path=path+(k,)))
                                        for k, v in val.items()])
    else:
        return val
# end: def transform_json_data(val, node_handler, path=())

def json_gather_leaf_jpaths(json_data):
    """Construct a map which for each leaf maps its jpath to the leaf value"""
    jpath2leaf = collections.OrderedDict()
    def save_leaf_path(val, path, is_leaf):
        if not is_leaf: 
            return val
        assert path not in jpath2leaf
        jpath2leaf[path] = val
        return val
    transform_json_data(json_data, save_leaf_path)
    return jpath2leaf

def map_vals(d):
    """Return an iterator over map values"""
    return map(operator.itemgetter(1), d.items())

def wraps(f):
    """Like functools.wraps but sets __wrapped__ even on Python 2.7"""
    wrapper = functools.wraps(f)
    wrapper.__wrapped__ = f
    return wrapper

def unwrap(f):
    """Find the original function under layers of wrappers"""
    return f if not hasattr(f, '__wrapped__') else unwrap(f.__wrapped__)

# UUID_RE: a regexp for UUIDs
def _uuid_part_re(length):
    return '[0-9a-f]{' + str(length) + '}'
UUID_RE = '-'.join(map(_uuid_part_re, (8, 4, 4, 4, 12)))

class MultilineFormatter(logging.Formatter):
    # adapted from https://stackoverflow.com/questions/45216826/how-to-print-multiline-logs-using-python-logging-module
    def format(self, record):
        """Break up long lines, for easier reading in Emacs"""
        result = super(MultilineFormatter, self).format(record)
        LINE_LEN = int(os.environ.get('VIRAL_NGS_LOG_LINE_LEN', '0'))
        if LINE_LEN and len(result) > LINE_LEN:
            result = '\n'.join([('  ' if i > 0 else '') + result[i:i+LINE_LEN] for i in range(0, len(result), LINE_LEN)])
        return result

class NestedParserAction(argparse.Action):

    """Parses a list of args using a specified parser.

    Code adapted from argparse._AppendAction
    """

    def __init__(self,
                 option_strings,
                 dest,
                 nargs=None,
                 const=None,
                 default=None,
                 type=None,
                 choices=None,
                 required=False,
                 help=None,
                 metavar=None,
                 nested_parser=None):
        if nested_parser is None:
            raise ValueError('nested_parser must be specified')
        if nargs == 0:
            raise ValueError('nargs for append actions must be > 0; if arg '
                             'strings are not supplying the value to append, '
                             'the append const action may be more appropriate')
        if const is not None:
            raise ValueError('const should not be specified')
        super(NestedParserAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=nargs,
            const=const,
            default=default,
            type=type,
            choices=choices,
            required=required,
            help=help,
            metavar=metavar)
        self.nested_parser = nested_parser

    def __call__(self, parser, namespace, values, option_string=None):
        items = copy.copy(self._argparse_ensure_value(namespace, self.dest, []))
        log.debug('NESTED PARSING %s', values)
        items.append(self.nested_parser.parse_args(values))
        setattr(namespace, self.dest, items)

    @staticmethod
    def _argparse_ensure_value(namespace, name, value):
        """Ensure that `namespace` has the attr `name`, with default value `value`,
        and return the attri's value."""
        if getattr(namespace, name, None) is None:
            setattr(namespace, name, value)
        return getattr(namespace, name)

# end: class _NestedParserAction(argparse.Action)

if hasattr(re, 'fullmatch'):
    re_fullmatch = re.fullmatch
else:
    def re_fullmatch(regex, string, flags=0):
        """Emulate python-3.4 re.fullmatch()."""
        return re.match("(?:" + regex + r")\Z", string, flags=flags)

  
def maybe_wait_for_result(arg, timeout=None):
    """If arg is a Future, wait for its result, else just return it."""
    if isinstance(arg, concurrent.futures.Future):
        log.info('waiting for result of future %s', arg)
    return arg.result(timeout=timeout) if isinstance(arg, concurrent.futures.Future) else arg

def maps(obj, *keys):
    """Return True if `obj` is a mapping that contains al `keys`."""
    return isinstance(obj, collections.Mapping) and all(k in obj for k in keys)

_str_type = basestring if hasattr(builtins, 'basestring') else str

def is_str(obj):
    """Test if obj is a string type, in a python2/3 compatible way.
    From https://stackoverflow.com/questions/4232111/stringtype-and-nonetype-in-python3-x
    """
    return isinstance(obj, _str_type)

def pretty_print_json(json_dict, sort_keys=True):
    """Return a pretty-printed version of a dict converted to json, as a string."""
    return json.dumps(json_dict, indent=4, separators=(',', ': '), sort_keys=sort_keys)

def write_json(fname, **json_dict):
    util.file.dump_file(fname=fname, value=pretty_print_json(json_dict))

def _load_dict_sorted(d):
    return collections.OrderedDict(sorted(d.items()))

def json_loads(s):
    """Load json from a string, canonicalizing order of fields within records."""
    return json.loads(s.strip(), object_hook=_load_dict_sorted, object_pairs_hook=collections.OrderedDict)

def json_loadf(fname):
    """Load json from a filename."""
    return json_loads(util.file.slurp_file(fname))

def _kill_proc_tree(pid, sig=signal.SIGTERM, include_parent=True,
                   timeout=None, on_terminate=None):
    """Kill a process tree (including grandchildren) with signal
    "sig" and return a (gone, still_alive) tuple.
    "on_terminate", if specified, is a callabck function which is
    called as soon as a child terminates.

    From https://psutil.readthedocs.io/en/latest/
    """
    if pid == os.getpid():
        raise RuntimeError("I refuse to kill myself")
    parent = psutil.Process(pid)
    children = parent.children(recursive=True)
    if include_parent:
        children.append(parent)
    for p in children:
        p.send_signal(sig)
    gone, alive = psutil.wait_procs(children, timeout=timeout,
                                    callback=on_terminate)
    return (gone, alive)

def kill_proc_tree(proc):
    """Kill a popen process tree"""
    gone, alive = _kill_proc_tree(proc.pid)
    time.sleep(1)
    # send SIGKILL
    for p in (alive or ()):
        p.kill()
    psutil.wait_procs(alive)

def find_free_port():
    """Find a free TCP port."""
    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        s.bind(('', 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]
