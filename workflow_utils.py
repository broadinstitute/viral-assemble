#!/usr/bin/env python

# * Preamble
# ** module comment

"""
Utilities for dealing with workflows, including cloud workflows.

Commands in this module help with the following:
   - jointly version code and data using git-annex ( https://git-annex.branchable.com/ ) and DataLad ( https://www.datalad.org/ )
   - keep track of data located in various places (dnanexus, AWS, GCS) in a unified and uniform way
   - keep track of analyses that have been run, storing each analysis' details in a single self-contained 'analysis directory'
     within a git/git-annex repository.
   - re-running an analysis, or a group of analyses, on any cloud (regardless of where they were originally run),
     while applying specified modifications to the parameters.  In particular, this enables curating sets of benchmark
     analyses, represented as a git subdir containing a group of analyses, and re-running them with a given parameter
     modification.

In combination, these utilities enable much of what DNAnexus enables, but:
   - without forcing all analyses and data to reside in one vendor
   - providing a uniform interface that works the same across vendors and clouds, as well as analyses run on
     UGER or Odyssey or locally
   - allowing the use of spot/preemptible instances (which DNAnexus disallows)
   - allowing easy archiving/unarchiving of data, using git-annex commands
   - using all of git's normal mechanisms, such as branches and submodules, to organize things
"""

__author__ = "ilya@broadinstitute.org"
__commands__ = []

# ** imports

# *** built-ins

import platform
assert platform.python_version().startswith('2.7')  # dxpy requirement; should drop and use command-line

import argparse
import logging
import json
import subprocess
import os
import os.path
import shutil
import glob
import collections
import time
import getpass
import uuid
import functools
import operator
import builtins
import datetime
import hashlib
import pipes
import shlex
import sys
import urllib
import SimpleHTTPServer
import SocketServer
import traceback

# *** 3rd-party
import dxpy
import dxpy.bindings.dxfile_functions
import boto3
import uritools
import pytz

# *** intra-module
import util.cmd
import util.file
import util.misc

_log = logging.getLogger(__name__)

#logging.basicConfig(format="%(asctime)s - %(module)s:%(lineno)d:%(funcName)s - %(levelname)s - %(message)s")
_log.setLevel(logging.DEBUG)


# * Utils
# ** Generic utils

_str_type = getattr(builtins, 'basestring', 'str')
_long_type = getattr(builtins, 'long', 'int')
_scalar_types = (type(None), _str_type, _long_type, int, bool, float)
if hasattr(builtins, 'unicode'):
    _scalar_types += (unicode,)

def _is_str(obj):
    """Test if obj is a string type, in a python2/3 compatible way.
    From https://stackoverflow.com/questions/4232111/stringtype-and-nonetype-in-python3-x
    """
    return isinstance(obj, _str_type)

def _is_mapping(obj):
    return isinstance(obj, collections.Mapping)

def _maps(obj, *keys):
    return _is_mapping(obj) and all(k in obj for k in keys)

def _is_scalar(val):
    isinstance(val, _scalar_types)

def _ord_dict(*args):
    return collections.OrderedDict(args)

def _dict_rename_key(d, old_key, new_key):
    print('d.keys=', d.keys())
    assert new_key not in d
    d[new_key] = d[old_key]
    del d[old_key]
    return d

def _pretty_print_json(json_dict, sort_keys=True):
    """Return a pretty-printed version of a dict converted to json, as a string."""
    return json.dumps(json_dict, indent=4, separators=(',', ': '), sort_keys=sort_keys)

def _write_json(fname, **json_dict):
    util.file.dump_file(fname=fname, value=_pretty_print_json(json_dict))

def _quote(s):
    return pipes.quote(s) if hasattr(pipes, 'quote') else shlex.quote(s)

def _make_cmd(cmd, *args):
    print('ARGS=', args)
    return ' '.join([cmd] + [_quote(str(arg)) for arg in args if arg not in (None, '')])

def _run(cmd, *args):
    cmd = _make_cmd(cmd, *args)
    print('running command: ', cmd, 'cwd=', os.getcwd())
    beg_time = time.time()
    subprocess.check_call(cmd, shell=True)
    print('command succeeded in {}s: {}'.format(time.time()-beg_time, cmd))

def _run_succeeds(cmd, *args):
    try:
        _run(cmd, *args)
        return True
    except subprocess.CalledProcessError:
        return False

def _run_get_output(cmd, *args):
    cmd = _make_cmd(cmd, *args)
    print('running command: ', cmd)
    beg_time = time.time()
    output = subprocess.check_output(cmd, shell=True)
    print('command succeeded in {}s: {}'.format(time.time()-beg_time, cmd))
    return output.strip()


def _json_loads(s):
    return json.loads(s.strip(), object_pairs_hook=collections.OrderedDict)

def _json_loadf(fname):
    return _json_loads(util.file.slurp_file(fname))

def _run_get_json(cmd, *args):
    return _json_loads(_run_get_output(cmd, *args))

def _transform_parsed_json(val, node_handler, path=()):
    """Transform a parsed json structure, by replacing nodes for which `node_handler` returns
    an object different from `val`, with that object.  If `node_handler` has a named arg `path`,
    the path from the root will be passed as that arg.
    """
    recurse = functools.partial(_transform_parsed_json, node_handler=node_handler)
    node_handler_named_args = util.misc.getnamedargs(node_handler)
    node_handler_args = {}
    if 'path' in node_handler_named_args:
        node_handler_args['path'] = path
    handled_val = node_handler(val, **node_handler_args)
    if handled_val is not val:
        print('resolved', val, 'to', handled_val)
        return handled_val
    if isinstance(val, list):
        return [recurse(val=v, path=path+(i,)) for i, v in enumerate(val)]
    if isinstance(val, collections.Mapping):
        return _ord_dict(*[(k, recurse(val=v, path=path+(k,)))
                           for k, v in val.items()])
    return val


def _json_to_org(val, org_file, depth=1, heading='root'):
    """Transform a parsed json structure to org.
    """
    with open(org_file, 'w') as out:
        def _recurse(val, heading, depth):
            def _header(s): out.write('*'*depth + ' ' + str(s) + '\n')
            def _line(s): out.write(' '*depth + str(s) + '\n')
            out.write('*'*depth + ' ' + heading)
            if isinstance(val, list):
                out.write(' - list of ' + str(len(val)) + '\n')
                if len(val):
                    for i, v in enumerate(val):
                        _recurse(v, heading=str(i), depth=depth+2)
            elif isinstance(val, collections.Mapping):
                out.write(' - map of ' + str(len(val)) + '\n')
                if len(val):
                    for k, v in val.items():
                        _recurse(v, heading='_'+k+'_', depth=depth+2)
            else:
                out.write(' - ' + str(val) + '\n')
        _recurse(val=val, heading=heading, depth=depth)

def workflow_utils_init():
    """Install the dependencies: cromwell and dxpy and git-annex."""
    _run('conda install cromwell dxpy git-annex')


# ** AnalysisDir

class AnalysisDir(object):

    """A directory representing an analysis.   Typically a directory in a git-annex repository,
    with analysis metadata in standard format, and analysis files as git-annex files under the 
    analysis directory.

    An analysis directory can represent an analysis in various states: being prepared, submitted,
    finished; succeeded or failed.

    """
    
    METADATA_VERSION = "1.0"

    def __init__(self, path):
        self.path = path

    def is_analysis_dir(self):
        return os.path.isfile(os.path.join(self.path, 'metadata.json'))

# ** DNAnexus-related utils

@util.misc.memoize
def _dx_describe(dxid):
    return _run_get_json('dx', 'describe', '--verbose', '--details', '--json', dxid)

DX_URI_PFX='dx://'

def _url_to_dxid(url):
    return os.path.splitext(url[len(DX_URI_PFX):].split('/')[0])[0]

def _standardize_dx_url(url):
    dxid = _url_to_dxid(url)
    dx_descr = _dx_describe(_url_to_dxid(url))
    return DX_URI_PFX + dxid + '/' + urllib.pathname2url(dx_descr['name'])

# ** git-annex-related utils

def _git_annex_lookupkey(f):
    return _run_get_output('git', 'annex', 'lookupkey', f)

def _git_annex_get_metadata(key, field):
    ga_mdata = _run_get_json('git', 'annex', 'metadata', '--json', '--key='+key)
    return ga_mdata.get('fields', {}).get(field, [])

def _git_annex_set_metadata(key, field, val):
    _run('git', 'annex', 'metadata', '--key='+key, '-s', field + '=' + val)

def _git_annex_get_url_key(url):
    """Return the git-annex key for a url"""
    print('get-url-key', url)
    assert uritools.isuri(url)
    url = _minimize_url(url)
    assert uritools.isuri(url)

    url_parts = uritools.urisplit(url)
    url_parts = url_parts._replace(path=uritools.uriencode(url_parts.path, safe='/'))
    url = uritools.uriunsplit(url_parts)

    assert uritools.isuri(url)
    return 'URL--' + url

    # old way, maybe less abstraction-breaking but slower -- let git-annex compute the key:
    # with util.file.tmp_dir(dir=os.getcwd()) as tmp_d:
    #     f = os.path.join(tmp_d, 'f')
    #     _run('git', 'annex', 'fromkey', '--force', url, f)
    #     key = _git_annex_lookupkey(f)
    #     shutil.rmtree(tmp_d, ignore_errors=True)
    #     return key


# * json_to_org


def json_to_org(json_fname, org_fname=None):
    """Convert json to org"""
    if not org_fname:
        org_fname = os.path.splitext(json_fname)[0]+'.org'
    _json_to_org(val=_json_loadf(json_fname), org_file=org_fname)

def parser_json_to_org(parser=argparse.ArgumentParser()):
    parser.add_argument('json_fname', help='json file to import')
    parser.add_argument('org_fname', help='org file to output', nargs='?')
    util.cmd.attach_main(parser, json_to_org, split_args=True)

__commands__.append(('json_to_org', parser_json_to_org))


# * submit_analysis_wdl
# ** get_workflow_inputs


#
# given a dx analysis, run it locally
#

def get_workflow_inputs(workflow_name, docker_img):
    """Run womtool to get the inputs of the wdl workflow"""
    _run('docker run --rm ' + docker_img + ' tar cf - source/pipes/WDL > wdl.tar')
    _run('tar xvf wdl.tar')
    for f in glob.glob('source/pipes/WDL/workflows/*.wdl'):
        shutil.copy(f, '.')
    for f in glob.glob('source/pipes/WDL/workflows/tasks/*.wdl'):
        shutil.copy(f, '.')
    shutil.rmtree('source')
    os.unlink('wdl.tar')
    return _run_get_json('womtool inputs ' + workflow_name + '.wdl')
    

# ** _get_dx_val
def _get_dx_val(val, dx_files_dir):
    """Resolve a dx value: if it is a scalar, just return that;
    if it is a dx file, fetch the file, cache it locally, and return the path to the file."""

    #print('parsing val: ', val)
    util.file.mkdir_p(dx_files_dir)
    if isinstance(val, list):
        return [_get_dx_val(val=v, dx_files_dir=dx_files_dir) for v in val]
    if isinstance(val, collections.Mapping) and '$dnanexus_link' in val:
        link = val['$dnanexus_link']
        if isinstance(link, collections.Mapping) and 'analysis' in link:
            print('link is ', link)
            return _get_dx_val(dxpy.DXAnalysis(link['analysis']).describe()['output'][link['stage']+'.'+link['field']],
                               dx_files_dir)
        elif (_is_str(link) and link.startswith('file-')) or (isinstance(link, collections.Mapping) and 'id' in link and _is_str(link['id']) and link['id'].startswith('file-')):
            if _is_str(link) and link.startswith('file-'):
                dxid = link
            else:
                dxid = link['id']
            descr = dxpy.describe(dxid)
            dx_file = os.path.join(dx_files_dir, dxid) + '-' + descr['name']
            file_size = int(descr['size'])

            # see if the file is cached in git-annex
            ga_mdata = _run_get_json('git annex metadata --json --key=WORM-s0-m0--dx-' + dxid)
            if ga_mdata['fields']:
                ga_key = ga_mdata['fields']['ga_key'][0]
                _run('git annex get --key', ga_key)
                _run('git annex fromkey', ga_key, dx_file)

            if not os.path.isfile(dx_file) or os.path.getsize(dx_file) != file_size:
                print('fetching', dxid, 'to', dx_file)
                # TODO: check that there is enough free space (with some to spare)
                if os.path.isfile(dx_file):
                    os.unlink(dx_file)
                fs_info = os.statvfs(dx_files_dir)
                assert file_size < fs_info.f_bsize * fs_info.f_bavail
                dxpy.bindings.dxfile_functions.download_dxfile(dxid=dxid, filename=dx_file+'.fetching', show_progress=True)
                assert os.path.getsize(dx_file+'.fetching') == file_size
                os.rename(dx_file+'.fetching', dx_file)
                print('fetched', dxid, 'to', dx_file)
                print('curdir is ', os.getcwd())
                _run('git annex add', dx_file)
                ga_key = _run_get_output('git annex lookupkey ' + dx_file).strip()
                _run('git annex metadata -s', 'dxid+=' + dxid, dx_file)
                # record a mapping from the dxid to the git-annex key
                _run('git annex metadata', '--key=WORM-s0-m0-dx-'+dxid, '-s', 'ga_key='+ga_key)
                # register a URL that can be used to re-fetch this file from DNAnexus;
                # the URL is only valid if the 'run_dx_url_server' command is running.
                #_run('git annex registerurl ' + ga_key + ' ' + ' http://localhost:8080/dx/' + dxid)
        else:
            raise RuntimeError('Cannot parse dx link {}'.format(link))
        return dx_file
    else:
        return val
# end: def _get_dx_val(val, dx_files_dir = 'input_files')

# ** proper implementation

class GitAnnexTool(object):

    def get_metadata(self, key):
        """Return metadata associated with key, as a tuple of values.  If no metadata is associated
        with key, return empty tuple."""
        key_metadata = _run_get_json("git annex metadata --json --key='{}'".format(key))
        return tuple(key_metadata.get('fields', {}).get(key, ()))

class GitAnnexRepo(object):

    """A local git-annex repository, together with its configured cloud remotes.
    """

    def from_key(ga_key, repo_rel_path):
        """Set the file at `repo_rel_path` to point to the git-annex key `ga_key`."""
        pass

    def record_file_uri(file_uri):
        """Ensures that the given URI is recorded in this repo.

        Returns:
           A git-annex key for the contents of file_uri.  
        """
        pass

    pass

class GitAnnexRemote(object):

    def has_key(key):
        """Test whether the remote contains a given key"""
        raise NotImplementedError()

    def put_key(key):
        """Ensure the remote has content with given key."""
        raise NotImplementedError()

# ** _resolve_file_values

class FileResolver(object):
    """Abstract base class for recognizing and resolving values that point to files."""

    def __init__(self, input_files_dir):
        self.input_files_dir = input_files_dir
    
    def is_file_value(self, val): 
        """Test whether `val` represents a file that this FileResolver can handle"""
        raise NotImplementedError()

    def create_git_annex_symlink(self, file_val, git_file_path):
        """Given a value representing a file, create a git-annex link to the file under the relative path `git_file_path`.
        Returns the git-annex key representing the contents of `file_val`.
        """
        raise NotImplementedError()

class DxFileResolver(FileResolver):

    """FileResolver for DNAnexus links"""

    def is_file_value(self, val):
        """Test whether `val` represents a DNAnexus file"""
        return _is_mapping(val) and '$dnanexus_link' in val

    def _get_dx_id(self, file_val):
        """Given a `val` that represents a DNAnexus file (possibly indirectly), get the dx ID (file-xxxxxxx) for the file."""
        link = file_val['$dnanexus_link']
        if isinstance(link, collections.Mapping) and 'analysis' in link:
            # a link to the output of a particular stage of a particular DNAnexus analysis
            return self._get_dx_file_id(dxpy.DXAnalysis(link['analysis']).describe()['output'][link['stage']+'.'+link['field']])
        elif (_is_str(link) and link.startswith('file-')) or \
             (isinstance(link, collections.Mapping) and 'id' in link and \
              _is_str(link['id']) and link['id'].startswith('file-')):
            if _is_str(link) and link.startswith('file-'):
                return link
            else:
                return link['id']

    def create_git_annex_symlink(self, file_val, git_file_path):
        """Given a value representing a file, create a git-annex link to the file under the relative path `git_file_path`.
        Returns the git-annex key representing the contents of `file_val`.
        """
        dxid = self._get_dx_id(file_val)


        dx_descr = dxpy.describe(dxid)
        dx_uri = 'dx://' + dxid + '/' + descr['name']

        _run("git annex addurl --fast '{}' --file '{}' ".format(dx_uri, git_file_path))
        ga_key = _run_get_output("git annex lookupkey '{}'".format(git_file_path))

        def get_md5e_key(uri_key):
            uri_key_mdata = _run_get_json("git annex metadata --json --key='{}'".format(uri_key))
            if uri_key_mdata['fields']:
                md5e_keys = [k for k in uri_key_mdata['fields'].get('altKeys', ()) if k.startswith('MD5E-')]
                if md5e_keys:
                    assert len(md5e_keys) == 1
                    return md5e_keys[0]
            return None
        md5e_key = get_md5e_key(ga_key)
        if md5e_key:
            _run('git rm ' + git_file_path)
            _run("git annex fromkey '{}' '{}'".format(md5e_key, git_file_path))
            ga_key = md5e_key
        return ga_key
    # end: def create_git_annex_symlink(self, file_val, git_file_path):
# end: class DxFileResolver(FileResolver)

def _resolve_file_values(input_name, val, input_files_dir):
    """For workflow inputs that point to a file, ensure that the file is present in the filesystem in which the 
    analysis will be run.

    Args:
      input_name: the name of an input to a workflow
      val: the value of the input.  May be a scalar, or a composite value such as a list or a map.
           May be a value representing a file.
    """

    #Resolve a dx value: if it is a scalar, just return that;
    #if it is a dx file, fetch the file, cache it locally, and return the path to the file.

    #print('parsing val: ', val)
    util.file.mkdir_p(dx_files_dir)
    if isinstance(val, collections.Mapping) and '$dnanexus_link' in val:
        link = val['$dnanexus_link']
        if isinstance(link, collections.Mapping) and 'analysis' in link:
            print('link is ', link)
            return _get_dx_val(dxpy.DXAnalysis(link['analysis']).describe()['output'][link['stage']+'.'+link['field']],
                               dx_files_dir)
        elif (_is_str(link) and link.startswith('file-')) or (isinstance(link, collections.Mapping) and 'id' in link and _is_str(link['id']) and link['id'].startswith('file-')):
            if _is_str(link) and link.startswith('file-'):
                dxid = link
            else:
                dxid = link['id']
            descr = dxpy.describe(dxid)
            dx_file = os.path.join(dx_files_dir, dxid) + '-' + descr['name']
            file_size = int(descr['size'])

            # see if the file is cached in git-annex
            ga_mdata = _run_get_json('git annex metadata --json --key=WORM-s0-m0--dx-' + dxid)
            if ga_mdata['fields']:
                ga_key = ga_mdata['fields']['ga_key'][0]
                _run('git annex get --key ' + ga_key)
                _run('git annex fromkey ' + ga_key + ' ' + dx_file)

            if not os.path.isfile(dx_file) or os.path.getsize(dx_file) != file_size:
                print('fetching', dxid, 'to', dx_file)
                # TODO: check that there is enough free space (with some to spare)
                if os.path.isfile(dx_file):
                    os.unlink(dx_file)
                fs_info = os.statvfs(dx_files_dir)
                assert file_size < fs_info.f_bsize * fs_info.f_bavail
                dxpy.bindings.dxfile_functions.download_dxfile(dxid=dxid, filename=dx_file+'.fetching', show_progress=True)
                assert os.path.getsize(dx_file+'.fetching') == file_size
                os.rename(dx_file+'.fetching', dx_file)
                print('fetched', dxid, 'to', dx_file)
                print('curdir is ', os.getcwd())
                _run('git annex add ' + dx_file)
                ga_key = _run_get_output('git annex lookupkey ' + dx_file).strip()
                _run('git annex metadata -s dxid+=' + dxid + ' ' + dx_file)
                # record a mapping from the dxid to the git-annex key
                _run('git annex metadata --key=WORM-s0-m0--{} -s ga_key={}'.format('dx-'+dxid, ga_key))
                # register a URL that can be used to re-fetch this file from DNAnexus;
                # the URL is only valid if the 'run_dx_url_server' command is running.
                _run('git annex registerurl ' + ga_key + ' ' + ' http://localhost:8080/dx/' + dxid)
        else:
            raise RuntimeError('Cannot parse dx link {}'.format(link))
        return dx_file
    else:
        return val
# end: def _resolve_file_values(val, dx_files_dir = 'input_files')


def _parse_cromwell_output_str(cromwell_output_str):
    """Parse cromwell output"""
    assert cromwell_output_str.count('Final Outputs:') == 1
    json_beg = cromwell_output_str.index('Final Outputs:') + len('Final Outputs:')
    json_end = cromwell_output_str.index('\n}\n', json_beg) + 2
    return _json_loads(cromwell_output_str[json_beg:json_end])

# ** import_dx_analysis


def _resolve_dx_link_to_dx_file_id_or_value(val, dx_analysis_id):
    """If `val` represents a DNAnexus link to a file, return {$dnanexus_link: file-xxxx}.
    If `val` represents a DNAnexus link to a value, return {$dnanexus_val: value}.
    Else, return None.
    """

    #print('parsing val: ', val)
    recurse = functools.partial(_resolve_dx_link_to_dx_file_id_or_value,
                                dx_analysis_id=dx_analysis_id)
    if not _maps(val, '$dnanexus_link'):
        return val

    link = val['$dnanexus_link']
    if _maps(link, 'stage') and ('field' in link or 'outputField' in link):
        print('link is ', link)
        linked_analysis_descr = _dx_describe(link.get('analysis', dx_analysis_id))
        linked_field = link['field'] if 'field' in link else link['outputField']
        return recurse(val=linked_analysis_descr['output'][link['stage']+'.'+linked_field])
    elif (_is_str(link) and link.startswith('file-')) or \
         (isinstance(link, collections.Mapping) and 'id' in link \
          and _is_str(link['id']) and link['id'].startswith('file-')):
        if _is_str(link) and link.startswith('file-'):
            dxid = link
        else:
            dxid = link['id']
        return {'$dnanexus_link': dxid}
    raise RuntimeError('Unknown $dnanexus_link: {}'.format(val))

def _resolve_link_dx(val, git_file_dir, dx_analysis_id):
    val_resolved = _resolve_dx_link_to_dx_file_id_or_value(val=val, dx_analysis_id=dx_analysis_id)
    if _maps(val_resolved, '$dnanexus_link'):
        dx_file_id = val_resolved['$dnanexus_link']
        assert dx_file_id.startswith('file-')
        util.file.mkdir_p(git_file_dir)
        val_resolved = {'$git_link': import_from_url(url='dx://' + dx_file_id, git_file_path=git_file_dir, fast=False)}
    return val_resolved

def _resolve_link(val, git_file_dir, methods):
    for method in methods:
        result = method(val=val, git_file_dir=git_file_dir)
        if result is not val: return result
    return val

def _resolve_links_in_parsed_json(val, rel_to_dir, methods, relpath='files'):
    """Given a parsed json structure, replace in it references to files (in various forms) with one uniform
    representation, a 'git link'.  A git link contains a relative path (relative to `analysis_dir`)
    pointing to a git file, typically a git-annex file.
    """
    def handle_node(val, path):
        maybe_resolve_link = _resolve_link(val=val,
                                           git_file_dir=os.path.join(rel_to_dir, relpath,
                                                                     *map(str, path)), methods=methods)
        if _maps(maybe_resolve_link, '$git_link'):
            maybe_resolve_link = {'$git_link': os.path.relpath(maybe_resolve_link['$git_link'], rel_to_dir)}
        return maybe_resolve_link
        
    return _transform_parsed_json(val, node_handler=handle_node)

def _resolve_links_in_dx_analysis(dx_analysis_id, analysis_dir):
    analysis_descr = _dx_describe(dx_analysis_id)
    del analysis_descr['originalInput']
    methods = [functools.partial(_resolve_link_dx, dx_analysis_id=dx_analysis_id)]
    resolved = _resolve_links_in_parsed_json(val=analysis_descr, rel_to_dir=analysis_dir, methods=methods)
    _write_json(os.path.join(analysis_dir, 'dx_resolved.json'), **resolved)

# def _resolve_link_dx_old(val, dx_analysis_id):
#     """Resolve a dx value: if it is a scalar, just return that;
#     if it is a dx link, resolve the link to a file-xxxx identifier."""

#     print('parsing val: ', val)
#     if not (_is_mapping(val) and '$dnanexus_link' in val):
#         return None

#     recurs = functools.partial(_resolve_dx_ids_in_val, analysis_descr=analysis_descr)
#     if isinstance(val, list):
#         return [recurs(val=v, git_path = os.path.join(git_path, str(i))) for i, v in enumerate(val)]
#     if isinstance(val, collections.Mapping):
#         if  '$dnanexus_link' not in val:
#             return collections.OrderedDict([(k, recurs(val=v, git_path=os.path.join(git_path, util.file.string_to_file_name(k))))
#                                             for k, v in val.items()])
#         link = val['$dnanexus_link']
#         if isinstance(link, collections.Mapping) and 'stage' in link and ('field' in link or 'outputField' in link):
#             print('link is ', link)
#             linked_analysis_descr = dxpy.DXAnalysis(link['analysis']).describe() if 'analysis' in link else analysis_descr
#             linked_field = link['field'] if 'field' in link else link['outputField']
#             return recurs(val=linked_analysis_descr['output'][link['stage']+'.'+linked_field],
#                           analysis_descr=analysis_descr, git_path=git_path)
#         elif (_is_str(link) and link.startswith('file-')) or (isinstance(link, collections.Mapping) and 'id' in link and _is_str(link['id']) and link['id'].startswith('file-')):
#             if _is_str(link) and link.startswith('file-'):
#                 dxid = link
#             else:
#                 dxid = link['id']




#     def _resolve_dx_ids_in_val(val, analysis_descr, git_path):
#         """Resolve a dx value: if it is a scalar, just return that;
#         if it is a dx link, resolve the link to a file-xxxx identifier."""

#         print('parsing val: ', val)

#         recurs = functools.partial(_resolve_dx_ids_in_val, analysis_descr=analysis_descr)
#         if isinstance(val, list):
#             return [recurs(val=v, git_path = os.path.join(git_path, str(i))) for i, v in enumerate(val)]
#         if isinstance(val, collections.Mapping):
#             if  '$dnanexus_link' not in val:
#                 return collections.OrderedDict([(k, recurs(val=v, git_path=os.path.join(git_path, util.file.string_to_file_name(k))))
#                                                 for k, v in val.items()])
#             link = val['$dnanexus_link']
#             if isinstance(link, collections.Mapping) and 'stage' in link and ('field' in link or 'outputField' in link):
#                 print('link is ', link)
#                 linked_analysis_descr = dxpy.DXAnalysis(link['analysis']).describe() if 'analysis' in link else analysis_descr
#                 linked_field = link['field'] if 'field' in link else link['outputField']
#                 return recurs(val=linked_analysis_descr['output'][link['stage']+'.'+linked_field],
#                               analysis_descr=analysis_descr, git_path=git_path)
#             elif (_is_str(link) and link.startswith('file-')) or (isinstance(link, collections.Mapping) and 'id' in link and _is_str(link['id']) and link['id'].startswith('file-')):
#                 if _is_str(link) and link.startswith('file-'):
#                     dxid = link
#                 else:
#                     dxid = link['id']
#                 util.file.mkdir_p(git_path)
#                 git_file_path = import_from_url(url='dx://' + dxid, git_file_path = git_path)
#                 return _ord_dict(('$git_relpath', os.path.relpath(git_file_path, )) # return { $git_link: }
#                 descr = dxpy.describe(dxid)
#                 dx_file = os.path.join(dx_files_dir, dxid) + '-' + descr['name']
#                 file_size = int(descr['size'])

#                 # see if the file is cached in git-annex
#                 ga_mdata = _run_get_json('git annex metadata --json --key=WORM-s0-m0--dx-' + dxid)
#                 if ga_mdata['fields']:
#                     ga_key = ga_mdata['fields']['ga_key'][0]
#                     _run('git annex get --key', ga_key)
#                     _run('git annex fromkey', ga_key, dx_file)

#                 if not os.path.isfile(dx_file) or os.path.getsize(dx_file) != file_size:
#                     print('fetching', dxid, 'to', dx_file)
#                     # TODO: check that there is enough free space (with some to spare)
#                     if os.path.isfile(dx_file):
#                         os.unlink(dx_file)
#                     fs_info = os.statvfs(dx_files_dir)
#                     assert file_size < fs_info.f_bsize * fs_info.f_bavail
#                     dxpy.bindings.dxfile_functions.download_dxfile(dxid=dxid, filename=dx_file+'.fetching', show_progress=True)
#                     assert os.path.getsize(dx_file+'.fetching') == file_size
#                     os.rename(dx_file+'.fetching', dx_file)
#                     print('fetched', dxid, 'to', dx_file)
#                     print('curdir is ', os.getcwd())
#                     _run('git annex add', dx_file)
#                     ga_key = _run_get_output('git annex lookupkey ' + dx_file).strip()
#                     _run('git annex metadata -s', 'dxid+=' + dxid, dx_file)
#                     # record a mapping from the dxid to the git-annex key
#                     _run('git annex metadata', '--key=WORM-s0-m0-dx-'+dxid, '-s', 'ga_key='+ga_key)
#                     # register a URL that can be used to re-fetch this file from DNAnexus;
#                     # the URL is only valid if the 'run_dx_url_server' command is running.
#                     #_run('git annex registerurl ' + ga_key + ' ' + ' http://localhost:8080/dx/' + dxid)
#             else:
#                 raise RuntimeError('Cannot parse dx link {}'.format(link))
#             return dx_file
#         else:
#             return val
# # end: def _resolve_dx_val(val, analysis_descr, dx_files_dir = 'input_files')




def _record_dx_metadata(val, analysis_dir, root_dir):
    """If 'val' is a filename, return a dict representing the file and some metadata about it;
    otherwise, return val as-is."""
    if isinstance(val, list): return [_record_file_metadata(v, analysis_dir, root_dir) for v in val]
    if _is_mapping(val): return collections.OrderedDict([(k, _record_file_metadata(v, analysis_dir, root_dir))
                                                          for k, v in val.items()])
    if not (_is_str(val) and (os.path.isfile(val) or os.path.isdir(val))): return val
    file_info = collections.OrderedDict([('_is_file' if os.path.isfile(val) else '_is_dir', True)])
    assert val.startswith(analysis_dir) or val.startswith(root_dir)
    if val.startswith(analysis_dir):
        relpath = os.path.relpath(val, analysis_dir)
        abspath = os.path.join(analysis_dir, relpath)
    else:
        cromwell_executions_dir = os.path.dirname(os.path.dirname(root_dir))
        relpath = os.path.relpath(val, cromwell_executions_dir)
        abspath = os.path.join(analysis_dir, 'output',
                               'call_logs' if os.path.basename(val) in ('stdout', 'stderr') else 'outputs', relpath)
        if os.path.isfile(val) and not os.path.isfile(abspath):
            print('LINKING {} to {}'.format(val, abspath))
            util.file.mkdir_p(os.path.dirname(abspath))
            shutil.copy(val, abspath)
        if os.path.isdir(val):
            util.file.mkdir_p(abspath)

    assert os.path.isabs(abspath) and abspath.startswith(analysis_dir), \
        'bad abspath: {} analysis_dir: {}'.format(abspath, analysis_dir)
    relpath = os.path.relpath(abspath, analysis_dir)
    assert not os.path.isabs(relpath), 'should be relative: {}'.format(relpath)
    assert os.path.isfile(abspath) or os.path.isdir(abspath), 'not file or dir: {}'.format(abspath)
    assert os.path.isdir(abspath) or os.path.getsize(abspath) == os.path.getsize(val)
    file_info['relpath'] = relpath
    if os.path.isfile(abspath):
        file_info['size'] = os.path.getsize(abspath)
        file_info['md5'] = _run_get_output('md5sum ' + abspath).strip().split()[0]
    return file_info


def import_dx_analysis(dx_analysis_id, analysis_dir_pfx):
    """Import a dnanexus analysis into git, in our format."""

    analysis_dir = analysis_dir_pfx + 'dx-' + dx_analysis_id
    util.file.mkdir_p(analysis_dir)
    mdata = _dx_describe(dx_analysis_id)

    methods = [functools.partial(_resolve_link_dx, dx_analysis_id=dx_analysis_id)]
    mdata = _resolve_links_in_parsed_json(val=mdata, rel_to_dir=analysis_dir, methods=methods)

    # convert to Cromwell's metadata format:
    # https://cromwell.readthedocs.io/en/develop/api/RESTAPI/#workflowmetadataresponse

    mdata['_metadata_version'] = AnalysisDir.METADATA_VERSION
    mdata['labels'] = _ord_dict(('platform', 'dnanexus'))

    _dict_rename_key(mdata, 'input', 'inputs')
    _dict_rename_key(mdata, 'output', 'outputs')
    _dict_rename_key(mdata, 'state', 'status')
    del mdata['originalInput']  # same as 'input'

    if mdata['status'] == 'done':
        mdata['status'] = 'Succeeded'
    elif mdata['status'] == 'failed':
        mdata['status'] = 'Failed'

    submission_datetime = datetime.datetime.fromtimestamp(float(mdata['created']) / 1000.0)
    tz_eastern = pytz.timezone('US/Eastern')
    mdata['submission'] = tz_eastern.localize(submission_datetime).isoformat('T')

    _write_json(os.path.join(analysis_dir, 'metadata.json'), **mdata)

def parser_import_dx_analysis(parser=argparse.ArgumentParser()):
    parser.add_argument('dx_analysis_id', help='dnanexus analysis id')
    parser.add_argument('analysis_dir_pfx', help='analysis dir prefix; analysis id will be added to it.')
    util.cmd.attach_main(parser, import_dx_analysis, split_args=True)

__commands__.append(('import_dx_analysis', parser_import_dx_analysis))

# ** submit_analysis_wdl

def submit_analysis_wdl(workflow_name, analysis_inputs_from_dx_analysis, docker_img, analysis_dir_pfx,
                        analysis_inputs_specified=None, 
                        analysis_descr='', analysis_labels=None,
                        data_repo=None, data_remote=None,
                        cromwell_server_url='http://localhost:8000',
                        backend='Local'):
    """Submit a WDL analysis to a Cromwell server.

    Args:
        workflow_name: name of the workflow, from pipes/WDL/workflows
        docker_img: docker image from which all viral-ngs and WDL code is taken.
        analysis_inputs_from_dx_analysis: id of a DNAnexus analysis from which to take analysis inputs
        analysis_inputs_specified: json file specifying analysis inputs directly; overrides any given by
         analysis_inputs_from_dx_analysis
        analysis_labels: json file specifying any analysis labels
        analysis_dir_pfx: prefix for the analysis dir
    """
    assert not os.path.isabs(analysis_dir)

    if not analysis_inputs_specified:
        analysis_inputs_specified = {}
    else:
        print('analysis_inputs_specified=', util.file.slurp_file(analysis_inputs_specified).strip())
        analysis_inputs_specified = _json_loadf(analysis_inputs_specified)

    analysis_id = create_analysis_id(workflow_name)
    print('ANALYSIS_ID is ', analysis_id)

    _run('docker pull ' + docker_img)
    docker_img_hash = docker_img + '@' + get_docker_hash(docker_img)

    if analysis_inputs_from_dx_analysis:
        dx_analysis = dxpy.DXAnalysis(dxid=analysis_inputs_from_dx_analysis)
        dx_analysis_descr = dx_analysis.describe()
    else:
        dx_analysis_descr = collections.OrderedDict(runInput={}, originalInput={})

    assert os.path.exists('.git/annex')
    data_repo = data_repo or os.getcwd()

    #_run('git config annex.security.allowed-http-addresses "127.0.0.1 ::1 localhost"')
                # if data_remote:
                #     _run('git annex enableremote ' + data_remote)

    t_dir = os.path.abspath(os.path.join(analysis_dir_pfx + '-' + analysis_id))
    util.file.mkdir_p(t_dir)
    with util.file.pushd_popd(t_dir):
        print('TTTTTTTTTTT t_dir=', t_dir)

        output_dir = os.path.join(t_dir, 'output')

        input_files_dir = os.path.join(t_dir, 'input_files')
        util.file.mkdir_p(input_files_dir)

        get_dx_val = functools.partial(_get_dx_val, dx_files_dir=input_files_dir)

        wdl_wf_inputs = get_workflow_inputs(workflow_name, docker_img=docker_img_hash)

        # TODO: use git annex batch mode to determine the keys for all the file-xxxx files, then
        # use batch mode to get the keys.  use -J to parallelize.  also option to use a cache remote, as
        # described at https://git-annex.branchable.com/tips/local_caching_of_annexed_files/

        # TODO: diff stages may have same-name inputs
        # note that dx runInput and originalInput mean the reverse of what one would think:
        # https://wiki.dnanexus.com/api-specification-v1.0.0/workflows-and-analyses
        dx_wf_inputs = { k.split('.')[-1] : v for k, v in dx_analysis_descr['runInput'].items()}   # check for dups
        dx_wf_orig_inputs = { k.split('.')[-1] : v for k, v in dx_analysis_descr['originalInput'].items()}   # check for dups
        print('DX_WF_RUNINPUTS', '\n'.join(dx_wf_inputs.keys()))
        print('DX_WF_ORIGINPUTS', '\n'.join(dx_wf_orig_inputs.keys()))
        new_wdl_wf_inputs = collections.OrderedDict()
        for wdl_wf_input, wdl_wf_input_descr in wdl_wf_inputs.items():
            wdl_wf_input_full = wdl_wf_input
            wdl_wf_input = wdl_wf_input.split('.')[-1]

            def _set_input(dx_wf_input):
                print('AAAAAAAAAAAAA wdl_wf_input=', wdl_wf_input, ' descr=', wdl_wf_input_descr, ' dx_wf_input=', dx_wf_input)
                new_wdl_wf_inputs[wdl_wf_input_full] = get_dx_val(dx_wf_input)

            if wdl_wf_input in analysis_inputs_specified:
                _set_input(analysis_inputs_specified[wdl_wf_input])
            elif wdl_wf_input in dx_wf_inputs:
                print('HAVE', wdl_wf_input, wdl_wf_input_descr, dx_wf_inputs[wdl_wf_input])
                _set_input(dx_wf_inputs[wdl_wf_input])
            elif (True or '(optional' not in wdl_wf_input_descr) and wdl_wf_input in dx_wf_orig_inputs:
                print('HAVE', wdl_wf_input, wdl_wf_input_descr, dx_wf_orig_inputs[wdl_wf_input])
                _set_input(dx_wf_orig_inputs[wdl_wf_input])
            else:
                print('MISSING', wdl_wf_input, wdl_wf_input_descr)
                assert '(optional' in wdl_wf_input_descr, \
                    'Missing required argument: {} {}'.format(wdl_wf_input, wdl_wf_input_descr)


        print(_pretty_print_json(new_wdl_wf_inputs))

        ################# put in the right docker ID!!  and find a place to keep the docker cache.

        def _rewrite_paths_to_gs(val):
            if _is_mapping(val): return {k: _rewrite_paths_to_gs(v) for k, v in val.items()}
            if isinstance(val, list): return list(map(_rewrite_paths_to_gs, val))
            if _is_str(val) and os.path.isfile(val):
                assert val.startswith(input_files_dir)
                val = 'gs://sabeti-ilya-cromwell/cromwell-inputs/' + analysis_id + '/' + val[len(input_files_dir)+1:]
            return val
        #_run('gsutil cp -r ' + input_files_dir + '/* gs://sabeti-ilya-cromwell/cromwell-inputs/' + analysis_id + '/')

        _write_json('inputs.json', **new_wdl_wf_inputs)
#            json.dump(_rewrite_paths_to_gs(new_wdl_wf_inputs), wf_out, indent=4, separators=(',', ': '))

        # TODO: option to update just some of the tasks.
        # actually, when compiling WDL, should have this option -- or, actually,
        # should make a new workflow where older apps are reused for stages that have not changed.
        _run('sed -i -- "s|{}|{}|g" *.wdl'.format('quay.io/broadinstitute/viral-ngs', docker_img_hash))

        util.file.mkdir_p(output_dir)
        util.file.mkdir_p(os.path.join(output_dir, 'outputs'))
        util.file.mkdir_p(os.path.join(output_dir, 'logs'))
        util.file.mkdir_p(os.path.join(output_dir, 'call_logs'))
        #util.file.mkdir_p(os.path.join(output_dir, 'metadata'))
        if backend == 'Local':
            wf_opts_dict = { "final_workflow_outputs_dir": os.path.join(output_dir, 'outputs'),
                             "final_workflow_log_dir": os.path.join(output_dir, 'logs'),
                             "final_call_logs_dir": os.path.join(output_dir, 'call_logs'),
                             "backend": "Local"
            }
        elif backend == 'JES':
            wf_opts_dict = {
                "final_workflow_log_dir": os.path.join(output_dir, 'logs'),
                "backend": "JES"
            }
        else:
            raise RuntimeError('Unknown backend - ' + backend)
        _write_json('cromwell_opts.json', **wf_opts_dict)
        _write_json('execution_env.json', ncpus=util.misc.available_cpu_count())

        _write_json('analysis_labels.json',
                    analysis_descr=analysis_descr,
                    docker_img=docker_img,
                    docker_img_hash=docker_img_hash,
                    inputs_from_dx_analysis=analysis_inputs_from_dx_analysis,
                    analysis_id=analysis_id,
                    analysis_dir=t_dir,
                    submitter=getpass.getuser(),
                    **dict(analysis_labels or {}))

        # add cromwell labels: dx project, the docker tag we ran on, etc.

        _log.info('Validating workflow')
        _run('womtool', 'validate',  '-i',  'inputs.json', workflow_name + '.wdl')
        _log.info('Validated workflow; calling cromwell')
        _run('zip imports.zip *.wdl')
        try:
            cromwell_output_str = _run_get_output('cromwell', 'submit', workflow_name+'.wdl',
                                                  '-t', 'wdl', '-i', 'inputs.json', '-l', 'analysis_labels.json',
                                                  '-o', 'cromwell_opts.json',
                                                  '-p', 'imports.zip', '-h', cromwell_server_url)
            cromwell_returncode = 0
        except subprocess.CalledProcessError as called_process_error:
            cromwell_output_str = called_process_error.output
            cromwell_returncode = called_process_error.returncode

        _log.info('Cromwell returned with return code %d', cromwell_returncode)
        util.file.dump_file('cromwell_submit_output.txt', cromwell_output_str)
        print('cromwell output is ', cromwell_output_str)
#        import sys
#        sys.exit(0)


        # util.file.dump_file(os.path.join(output_dir, 'cromwell_output.txt'), cromwell_output_str)
        # _run('sed -i -- "s|{}|{}|g" cromwell_execution_metadata.json'.format(t_dir+'/', ''))

        # if cromwell_returncode == 0:
        #     def make_paths_relative(v):
        #         if _is_str(v) and os.path.isabs(v) and v.startswith(t_dir):
        #             return os.path.relpath(v, t_dir)
        #         if isinstance(v, list):
        #             return list(map(make_paths_relative, v))
        #         return v
        #     cromwell_output_json = {k: make_paths_relative(v)
        #                             for k, v in _parse_cromwell_output_str(cromwell_output_str).items()}
        #     util.file.dump_file('outputs.json',
        #                         _pretty_print_json(cromwell_output_json))
        #     util.file.make_empty('analysis_succeeded.txt')
        # else:
        #     util.file.make_empty('analysis_failed.txt')

        _run('rm imports.zip *.wdl')
        if cromwell_returncode:
            raise RuntimeError('Cromwell failed - ' + cromwell_output_str)

#                 _run('sudo chown -R $USER . || true')
#                 _run('git annex add')
#                 _run('git commit -m after_running_analysis_' + analysis_id + '.')
# #                _run('git annex initremote content type=directory directory=/ndata/git-annex-content/ encryption=none')
#                 _run('git annex move --all --to {} -J{}'.format(data_remote or 'origin', util.misc.available_cpu_count()))
#                 _run('git annex dead here')
#                 # pull and merge here first? and try the sync several times, after a pause maybe
#                 _run('git annex sync --message git_annex_sync_analysis_{}'.format(analysis_id))

#                 # enable cleanup
#                 _run('chmod -R u+w . || true')


def create_analysis_id(workflow_name):
    """Generate a unique ID for the analysis."""
    return util.file.string_to_file_name('-'.join(map(str, 
                                                       ('analysis', time.strftime('%Y%m%d-%H%M%S', time.localtime())[2:], 
                                                        uuid.uuid4(), workflow_name))))[:1024]

def get_docker_hash(docker_img):
    if docker_img.startswith('sha256:'):
        return docker_img
    digest_lines = _run_get_output("docker images --digests --no-trunc --format "
                                   "'{{.Repository}}:{{.Tag}} {{.Digest}}'")
    digest_lines = [line for line in digest_lines.strip().split('\n') if docker_img in line]
    assert len(digest_lines) == 1
    digest_line = digest_lines[0]
    print('digest_line is |||{}|||'.format(digest_line))
    img, digest = digest_line.strip().split()
    assert img == docker_img + (':latest' if ':' not in docker_img else '')
    assert digest.startswith('sha256:') and len(digest) == 71
    return digest

def hash_and_store_file(fname):
    '''Store the given file in a key-value store, '''
    pass

def gather_run_results(docker_img, cromwell_output):
    pass
    
    #asdf

    # so, we also need to generate a unique id for this run;
    # since there can be multiple runs.

    # also, need to see how to handle failures, full or partial.
    #
#     sdb_client = boto3.client('sdb')
#    response = sdb_client.put_attributes(DomainName='viral_ngs_benchmarks',
#                                         ItemName=create_run_id(),
#                                         Attributes=[
                                            
#                                         ])

# =======================

########################################################################################################################


def parser_submit_analysis_wdl(parser=argparse.ArgumentParser()):
    parser.add_argument('workflow_name', help='Workflow name')
    parser.add_argument('--analysisDirPfx', dest='analysis_dir_pfx', default='pipelines/an',
                        help='directory where analysis will be stored; a unique suffix will be added')
    parser.add_argument('--analysisInputsFromDxAnalysis', dest='analysis_inputs_from_dx_analysis',
                        help='DNAnexus analysis ID to take analysis inputs from; specific ones can be overridden by '
                        '--analysisInputsSpecified')
    parser.add_argument('--dockerImg', dest='docker_img', default='quay.io/broadinstitute/viral-ngs')
    parser.add_argument('--analysisInputsSpecified', dest='analysis_inputs_specified',
                        help='explicitly specified analysis inputs')
    parser.add_argument('--analysisDescr', dest='analysis_descr', default='an analysis', help='description of the run')
    parser.add_argument('--dataRepo', dest='data_repo', help='git data repository')
    parser.add_argument('--dataRemote', dest='data_remote', help='git-annex data remote')
    parser.add_argument('--analysisLabels', dest='analysis_labels', nargs=2, action='append',
                        help='labels to attach to the analysis')
    parser.add_argument('--backend', default='Locall', help='backend on which to run')
    util.cmd.attach_main(parser, submit_analysis_wdl, split_args=True)

__commands__.append(('submit_analysis_wdl', parser_submit_analysis_wdl))

########################################################################################################################

    

########################################################################################################################

# * RedirectToCloudObject

class RedirectToCloudObject(SimpleHTTPServer.SimpleHTTPRequestHandler):

   def _get_dx(self):
       try:
           file_info = _run_get_json('dx', 'describe', '--json', self._get_dx_id())
       except subprocess.CalledProcessError as e:
           print('CalledProcessError:', e)
           raise
       print('got file_info', file_info)
       assert file_info['id'] == self._get_dx_id()

       self.send_response(307)
       if False and 'media' in file_info:
           self.send_header('Content-Type', file_info['media'])
       if False and 'size' in file_info:
           self.send_header('Content-Length', str(file_info['size']))
       #self.send_header('Content-Disposition', 'attachment; filename="{}"'.format(file_info['name']))

       new_path = _run_get_output('dx make_download_url --duration 2h ' + self._get_dx_id())
       print('new_path=', new_path)
       self.send_header('Location', new_path)
       self.end_headers()


   def _get_gs(self):
       _log.info('IN_gs_gs')
       try:
           _log.info('gs_url is %s', self._get_gs_uri())
           signed_url = _run_get_output('gsutil signurl -d 10m ' + self.server.gs_key + ' ' + self._get_gs_uri()).strip().split()[-1]
           _log.info('signed url is %s', signed_url)
       except subprocess.CalledProcessError as e:
           _log.error('CalledProcessError: %s', e)
           raise

       self.send_response(307)
       if False and 'media' in file_info:
           self.send_header('Content-Type', file_info['media'])
       if False and 'size' in file_info:
           self.send_header('Content-Length', str(file_info['size']))
       #self.send_header('Content-Disposition', 'attachment; filename="{}"'.format(file_info['name']))

       self.send_header('Location', signed_url)
       self.end_headers()

   def do_GET(self):
       try:
           _log.info('path is %s', self.path)
           _log.info('is gs path? %s', self.path.startswith('/gs/'))
           if not any(self.path.startswith('/' + p) for p in self.server.valid_paths):
               raise IOError('not a valid path')
           if self.path.startswith('/dx/'):
               self._get_dx()
           elif self.path.startswith('/gs/') and self.server.gs_key:
               _log.info('CALLING get_gs')
               self._get_gs()
               _log.info('RETURNED FROM get_gs')
           else:
               print('UNKNOWN PATH!!', self.path)
               raise IOError('Unknown URL')
       except (IOError, subprocess.CalledProcessError):
           self.send_error(404, 'file not found')

   def do_HEAD(self):
       self.do_GET()

   def _get_dx_id(self):
       return self.path[len('/dx/'):]

   def _get_gs_uri(self):
       return 'gs://' + self.path[len('/gs/'):]

def run_cloud_object_url_server(port, gs_key, valid_paths):
    """Start a webserver that will redirect URL requests of the form http://localhost/dx/file-xxxxxx to DNAnexus,
    requests of the form http://localhost/gs/gsbucket/gsobjectpath to Google Cloud Storage, and requests
    of the form http://localhost/s3/s3bucket/gsobjectpath to AWS S3.
    This gives each cloud object a stable http URL, permitting such files to be added to git-annex with the addurl command.
    """

    assert valid_paths, 'Some valid paths must be specified'
    server = None
    try:           
        SocketServer.TCPServer.allow_reuse_address = True
        server = SocketServer.TCPServer(("", port), RedirectToCloudObject)
        server.gs_key = gs_key
        server.valid_paths = valid_paths
        server.serve_forever()
    except Exception as e:
        if server is None:
            print('No server!')
        else:
            print('calling shutdown...')
            server.shutdown()
            print('calling server_close...')
            server.server_close()
            print('re-raising exception', e)
            raise

def parser_run_cloud_object_url_server(parser=argparse.ArgumentParser()):
    parser.add_argument('--port', default=8080, help='Port on which to run the webserver')
    parser.add_argument('--gsKey', dest='gs_key', help='Key for signing gs urls')
    parser.add_argument('--validPaths', dest='valid_paths', required=True, help='paths that may be accessed',
                        nargs='+')
    util.cmd.attach_main(parser, run_cloud_object_url_server, split_args=True)

__commands__.append(('run_cloud_object_url_server', parser_run_cloud_object_url_server))

########################################################################################################################

# * finalize_analysis_dir

# ** CromwelServer

class CromwellServer(object):

    """Interactions with a running Cromwell server"""

    def __init__(self, host):
        self.host = host

    def _api(self, url):
        return _run_get_json('curl -s -X GET "http://{}/api/workflows/v1/{}" -H "accept: application/json"'.format(self.host, url))

    def get_workflows(self):
        return self._api('query')['results']
    
    def get_metadata(self, workflow_id):
        return self._api('{}/metadata?expandSubWorkflows=false'.format(workflow_id))

# end: class CromwellServer(object)

def _is_analysis_done(analysis_dir):
    return os.path.exists(os.path.join(analysis_dir, 'output', 'logs'))

def _record_file_metadata(val, analysis_dir, root_dir):
    """If 'val' is a filename, return a dict representing the file and some metadata about it;
    otherwise, return val as-is."""
    if isinstance(val, list): return [_record_file_metadata(v, analysis_dir, root_dir) for v in val]
    if _is_mapping(val): return collections.OrderedDict([(k, _record_file_metadata(v, analysis_dir, root_dir))
                                                          for k, v in val.items()])
    if not (_is_str(val) and (os.path.isfile(val) or os.path.isdir(val))): return val
    file_info = collections.OrderedDict([('_is_file' if os.path.isfile(val) else '_is_dir', True)])
    assert val.startswith(analysis_dir) or val.startswith(root_dir)
    if val.startswith(analysis_dir):
        relpath = os.path.relpath(val, analysis_dir)
        abspath = os.path.join(analysis_dir, relpath)
    else:
        cromwell_executions_dir = os.path.dirname(os.path.dirname(root_dir))
        relpath = os.path.relpath(val, cromwell_executions_dir)
        abspath = os.path.join(analysis_dir, 'output',
                               'call_logs' if os.path.basename(val) in ('stdout', 'stderr') else 'outputs', relpath)
        if os.path.isfile(val) and not os.path.isfile(abspath):
            print('LINKING {} to {}'.format(val, abspath))
            util.file.mkdir_p(os.path.dirname(abspath))
            shutil.copy(val, abspath)
        if os.path.isdir(val):
            util.file.mkdir_p(abspath)

    assert os.path.isabs(abspath) and abspath.startswith(analysis_dir), \
        'bad abspath: {} analysis_dir: {}'.format(abspath, analysis_dir)
    relpath = os.path.relpath(abspath, analysis_dir)
    assert not os.path.isabs(relpath), 'should be relative: {}'.format(relpath)
    assert os.path.isfile(abspath) or os.path.isdir(abspath), 'not file or dir: {}'.format(abspath)
    assert os.path.isdir(abspath) or os.path.getsize(abspath) == os.path.getsize(val)
    file_info['relpath'] = relpath
    if os.path.isfile(abspath):
        file_info['size'] = os.path.getsize(abspath)
        file_info['md5'] = _run_get_output('md5sum ' + abspath).strip().split()[0]
    return file_info

def is_analysis_dir(d):
    """Test whether a given directory is an analysis dir"""
    return all(os.path.isfile(os.path.join(d, f)) for f in ('analysis_labels.json', 'inputs.json', 'cromwell_opts.json'))

def _gather_analysis_dirs(analysis_dirs_roots, processing_stats):
    def _get_analysis_dirs(d):
        if not os.path.isdir(d): return []
        if not os.access(d, os.R_OK | os.X_OK):
            processing_stats['dir_skipped_no_access'] += 1
            return []
        if is_analysis_dir(d):
            return [d]
        return functools.reduce(operator.concat, [_get_analysis_dirs(os.path.join(d, subd))
                                                  for subd in os.listdir(d)], [])
    return functools.reduce(operator.concat, map(_get_analysis_dirs, util.misc.make_seq(analysis_dirs_roots)), [])

# ** finalize_analysis_dirs impl
def finalize_analysis_dirs(cromwell_host, analysis_dirs_roots):
    """After a submitted cromwell analysis has finished, save results to the analysis dir.
    Save metadata, mark final workflow result, make paths relative to analysis dir."""
    cromwell_server = CromwellServer(host=cromwell_host)
    processing_stats = collections.Counter()
    _log.info('Gathering analysis dirs...')
    all_analysis_dirs = _gather_analysis_dirs(analysis_dirs_roots, processing_stats)
    processing_stats['all_analysis_dirs'] = len(all_analysis_dirs)
    _log.info('Got %d analysis dirs', len(all_analysis_dirs))
    for wf in cromwell_server.get_workflows():
        processing_stats['workflowsFromCromwell'] += 1
        mdata = cromwell_server.get_metadata(wf['id'])
        assert mdata['id'] == wf['id']
        assert 'workflowLog' not in mdata or mdata['workflowLog'].endswith('workflow.{}.log'.format(wf['id']))
        if 'analysis_dir' in mdata['labels']:
            analysis_dir = mdata['labels']['analysis_dir']
        else:
            processing_stats['noAnalysisDirForWorkflow'] += 1
            continue

        mdata_fname = os.path.join(analysis_dir, 'metadata_orig.json') # mdata['workflowLog'][:-4]+'.metadata.json'
        mdata_rel_fname = os.path.join(analysis_dir, 'metadata.json') # mdata['workflowLog'][:-4]+'.metadata.json'
        if os.path.exists(mdata_fname):
            processing_stats['metadata_already_saved'] += 1
        elif 'workflowRoot' not in mdata:
            processing_stats['workflow_root_not_in_mdata'] += 1
        else:
            _write_json(mdata_fname, **mdata)
            mdata_rel = _record_file_metadata(mdata, analysis_dir, mdata['workflowRoot'])
            _write_json(mdata_rel_fname, **mdata_rel)
            _log.info('Wrote metadata to %s and %s', mdata_fname, mdata_rel_fname)
            processing_stats['saved_metata'] += 1
    _log.info('Processing stats: %s', str(processing_stats))

def parser_finalize_analysis_dirs(parser=argparse.ArgumentParser()):
    parser.add_argument('analysis_dirs_roots', nargs='+', help='where are analysis dirs found')
    parser.add_argument('--cromwellHost', dest='cromwell_host', default='localhost:8000', help='cromwell server hostname')
    util.cmd.attach_main(parser, finalize_analysis_dirs, split_args=True)

__commands__.append(('finalize_analysis_dirs', parser_finalize_analysis_dirs))

########################################################################################################################

# * import_from_url


def _standardize_url(url):
    if url.startswith(DX_URI_PFX):
        return _standardize_dx_url(url)
    return url

def _minimize_url(url):
    if url.startswith(DX_URI_PFX):
        url = DX_URI_PFX+_url_to_dxid(url)
    return url


_url2md5key = {}

def _get_url_md5_key(url_key):
    if url_key not in _url2md5key:
        alt_keys = [k for k in _git_annex_get_metadata(key=url_key, field='alt_keys')
                    if k.startswith('MD5')]
        if alt_keys:
            _url2md5key[url_key] = alt_keys[0]
    return _url2md5key.get(url_key, None)

def _set_url_md5_key(url_key, md5_key):
    if url_key in _url2md5key:
        assert _url2md5key[url_key] == md5_key
        return
    _url2md5key[url_key] = md5_key
    _git_annex_set_metadata(key=url_key, field='alt_keys', val=md5_key)
    

def import_from_url(url, git_file_path, fast=False):
    """Imports a URL into the git-annex repository."""

    url = _standardize_url(url)
    if os.path.isdir(git_file_path):
        git_file_path = os.path.join(git_file_path, os.path.basename(url))

    assert not os.path.exists(git_file_path)

    # check if we have an md5 key for the url.
    url_key = _git_annex_get_url_key(url)
    md5_key = _get_url_md5_key(url_key)
    if md5_key:
        _run('git', 'annex', 'fromkey', md5_key, git_file_path)
        if not fast:
            _run('git', 'annex', 'get', git_file_path)
        return git_file_path

    _run('git', 'annex', 'addurl', '--fast' if fast else None, '--file', git_file_path, url)
    if not fast:
        _set_url_md5_key(url_key=url_key, md5_key=_git_annex_lookupkey(git_file_path))

    return git_file_path
    # check that the state of a dx object is closed.
    # (have option to) still compute md5.  if we don't write to disk, it'll be fast?

def parser_import_from_url(parser=argparse.ArgumentParser()):
    parser.add_argument('url', help='the URL of the file to add.')
    parser.add_argument('--gitFilePath', dest='git_file_path', help='filename in the git-annex repository', default='.')
    parser.add_argument('--fast', default=False, action='store_true', help='do not immediately download the file')
    util.cmd.attach_main(parser, import_from_url, split_args=True)

__commands__.append(('import_from_url', parser_import_from_url))

########################################################################################################################

def _gs_stat(gs_url):
    return _run_succeeds('gsutil stat', gs_url)

def _git_annex_checkpresentkey(key, remote):
    return _run_succeeds('git annex checkpresentkey', key, remote)

def _md5_base64(s):
    return hashlib.md5(s).digest().encode('base64').strip()

def _get_gs_remote_uuid():
    return '0b42380a-45f8-4b9d-82b3-e10aaf7bab6c'  # TODO determine this from git and git-annex config

def _copy_to_gs(git_file_path, gs_prefix = 'gs://sabeti-ilya-cromwell'):   # TODO make gs prefix configurable
    """Ensure the given file exists in gs."""

    if gs_prefix.endswith('/'):
        gs_prefix = gs_prefix[:-1]
    whereis = _run_get_json('git annex whereis --json ' + git_file_path)
    locs = []
    assert whereis['success']
    urls = [url for w in whereis['whereis'] for url in w.get('urls', ()) if url.startswith(gs_prefix) and _gs_stat(url)]
    fname = os.path.basename(git_file_path)
    urls_with_right_fname = [url for url in urls if url.endswith(fname)]

    key = _run_get_output('git annex lookupkey', git_file_path)
    gs_remote_uuid = _get_gs_remote_uuid()
    if not urls_with_right_fname:
        gs_url = '/'.join((gs_prefix, 'inp', _md5_base64(whereis['key']), fname))
        _run('git annex registerurl', key, gs_url)
        headers = None
        # TODO: set Content-Type based on extension
        # TODO: maybe use -z to compress if large text file.

        # TODO: if have a dx URL, use the Storage Transfer Service to download directly from 
        # dx and avoid transfer charges, if we're not running on google here.
        # or, run a wdl workflow with the temp url as parameter.

        if urls:
            _run('gsutil cp', urls[0], gs_url)
        else:
            _run('git annex get', git_file_path)
            _run('gsutil cp', git_file_path, gs_url)
            _run('git annex setpresentkey', key, gs_remote_uuid, 1)
        urls_with_right_fname.append(gs_url)
    assert _gs_stat(urls_with_right_fname[0])
    assert _git_annex_checkpresentkey(key, gs_remote_uuid)
    return urls_with_right_fname[0]

def _stage_file_for_backend(git_file_path, backend):
    """Ensure file exists in given backend"""
    if backend == 'Local':
        _run('git annex get ' + git_file_path)
        return git_file_path
    elif backend == 'JES':
        return _copy_to_gs(git_file_path)



########################################################################################################################

def _load_workflow_metadata(analysis_dir):
    return _json_loadf(glob.glob(os.path.join(analysis_dir, 'output', 'logs', 'workflow.*.json'))[0])

def analysis_dir_to_json(analysis_dir):
    """Gather data from one analysis dir into simple json"""
    mdata = _load_workflow_metadata(analysis_dir)

def gather_analyses(dirs):
    """Gather analyses from various sources into one simple uniform format."""
    pass

########################################################################################################################


# * compare_analysis_pairs

def _load_analysis_metadata(analysis_dir):
    """Read the analysis metadata from an analysis dir"""
    md = _json_loadf(os.path.join(analysis_dir, 'metadata.json'))
    md['analysis_dir'] = analysis_dir
    return md

def _flatten(val, pfx=''):
    if isinstance(val, list): return functools.reduce(operator.concat,
                                                      [_flatten(v, pfx+('.' if pfx else '')+str(i)) for i, v in enumerate(val)],
                                                      [])
    if _is_mapping(val):
        if '_is_file' in val:
            if 'md5' in val: return _flatten('md5:'+val['md5'], pfx)
            # for now, assume that if the file is from a URL then its unique id e.g. dx file-xxxx is in filename
            return _flatten(val['relpath'], pfx)
        if '_is_dir' in val:
            return _flatten(val['relpath'], pfx)
        return functools.reduce(operator.concat,
                                [_flatten(v, pfx+('.' if pfx else '')+k) for k, v in val.items()],
                                [])
    return [(pfx, val)]

def compare_analysis_pairs(workflow_name, analysis_dirs, filter_A, filter_B, metrics):
    """Compare pairs of analyses from `analysis_dirs` where the two analyses of a pair have
    identical values of `common_fields`, the first analysis of the pair matches the criteria in `filter_A` and
    the second matches the criteria in `filter_B`.  For such pairs, show the distribution of output 
    `metrics`.
    """

    # For each distinct combination of `common_fields` values, gather the analyses with that combination.
    analysis_dirs=[os.path.join(analysis_dirs, d) for d in os.listdir(analysis_dirs) \
                   if os.path.isfile(os.path.join(analysis_dirs, d, 'metadata.json'))]
    filter_A = dict(filter_A)
    filter_B = dict(filter_B)
    print('filter_A=', filter_A)
    print('filter_B=', filter_B)
    metrics = metrics or ()

    flat_analyses = [ _flatten(_load_analysis_metadata(analysis_dir)) for analysis_dir in analysis_dirs if
                      os.path.isfile(os.path.join(analysis_dir, 'metadata.json')) ]
    flat_analyses = [ f for f in flat_analyses if dict(f)['status'] == 'Succeeded' ]
    print('len(analysis_dirs)=', len(analysis_dirs))
    print('len(flat_analyses)=', len(flat_analyses))

    assert filter_A.keys() == filter_B.keys(), 'Atypical usage -- probably error?'
    
    filter_fields = set(set(filter_A.keys()) | set(filter_B.keys()))

    def _get_common_fields(flat_analysis):
        return tuple(sorted([(k, v) for k, v in flat_analysis \
                             if k not in filter_fields and \
                             (k.startswith('inputs.') or k.startswith('labels.docker_img_hash'))]))

    analysis_groups = collections.defaultdict(list)
    for a in flat_analyses:
        analysis_groups[_get_common_fields(a)].append(a)

    def _matches(a, filt):
        return all([k not in filt or filt[k] == v for k, v in a])

    metric2diffs = collections.defaultdict(list)

    found_pairs = 0
    for ag, analyses in analysis_groups.items():
        ag_A = [a for a in analyses if _matches(a, filter_A)]
        ag_B = [a for a in analyses if _matches(a, filter_B)]
        if ag_A and ag_B:
            an_A = collections.OrderedDict(ag_A[0])
            an_B = collections.OrderedDict(ag_B[0])
            #print('\nag=', ag, '\nan_A=', an_A, '\nan_B=', an_B)
            found_pairs += 1
            for metric in metrics:
                metric2diffs[metric].append((an_B[metric] - an_A[metric], an_A['analysis_dir'], an_B['analysis_dir']))
    print('found_pairs=', found_pairs)
    for metric in metrics:
        print('metric=', metric, 'diffs=\n')
        print('\n'.join(map(str, sorted(metric2diffs[metric]))))

def parser_compare_analysis_pairs(parser=argparse.ArgumentParser()):
    parser.add_argument('--workflowName', dest='workflow_name', required=True)
    parser.add_argument('--filterA', dest='filter_A', nargs=2, action='append', required=True)
    parser.add_argument('--filterB', dest='filter_B', nargs=2, action='append', required=True)
    parser.add_argument('--metrics', action='append')
    parser.add_argument('--analysisDirs', dest='analysis_dirs', default='pipelines', help='dir containing analysis dirs')
    util.cmd.attach_main(parser, compare_analysis_pairs, split_args=True)

__commands__.append(('compare_analysis_pairs', parser_compare_analysis_pairs))


# ** analyze_workflows_orig

def analyze_workflows_orig():
    bam2analyses = collections.defaultdict(list)
    stats = collections.Counter()
    for analysis in os.listdir('runs'):
        def _file(name): return os.path.join('runs', analysis, name)
        if not os.path.isfile(_file('cromwell_execution_metadata.json')):
            stats['no_metadata'] += 1
            continue
        mdata = _json_loadf(_file('cromwell_execution_metadata.json'))
        if not mdata['workflowName'].startswith('assemble_denovo'):
            stats['wrong_workflowName'] += 1
            continue
        if not mdata['status'] == 'Succeeded':
            stats['wrong_status'] += 1
            continue
        stats['considered'] += 1
        mdata = copy.deepcopy(mdata)
        mdata['analysis_id'] = analysis
        bam = mdata['inputs']['assemble_denovo.reads_unmapped_bam']
        bam2analyses[bam].append(mdata)
        print(analysis)
    print(stats)
    diffs = []
    for bam, mdatas in bam2analyses.items():
        if len(mdatas) > 1:
            print('bam=', bam, 'analyses=', '\n'.join(map(str, [mdata['analysis_id'] for mdata in mdatas])))

            def canonicalize_val(v, analysis_id):
                if not _is_str(v): return v
                orig_v = v
                v = os.path.join('runs', analysis_id, v)
                if _is_str(v) and os.path.islink(v) and os.path.lexists(v) and '.git/annex/objects' in os.path.realpath(v):
                    return os.path.basename(os.path.realpath(v))
                #print('not canon: ', v)
                return orig_v

            for i in range(len(mdatas)):
                for j in range(i+1, len(mdatas)):
                    if mdatas[i]['outputs']['assemble_denovo.refine_2x_and_plot.assembly_length_unambiguous'] == \
                       mdatas[j]['outputs']['assemble_denovo.refine_2x_and_plot.assembly_length_unambiguous']: continue
                    print('COMPARING ', mdatas[i]['analysis_id'], mdatas[j]['analysis_id'])
                    def show_diff(which):
                        def get_items(mdata):
                            analysis_id = mdata['analysis_id']
                            mdata = mdata[which]
                            mdata = {k: canonicalize_val(v, analysis_id) for k, v in mdata.items()}
                            mdata = {k:v for k,v in mdata.items() if isinstance(v, (int, float)) or _is_str(v) and not v.startswith('SHA256')}
                            #print('canonicalized version', mdata)
                            return mdata
                        items_i = get_items(mdatas[i])
                        items_j = get_items(mdatas[j])
                        print('\n\n\n')
                        #print('KEYS', items_i.keys())
                        assert items_i.keys() == items_j.keys()
                        print('------------differing items {}----------', which)
                        for k in sorted(items_i.keys()):
                            if items_i[k] != items_j[k]:
                                print(k, items_i[k], items_j[k])
                                if k == 'assemble_denovo.refine_2x_and_plot.assembly_length_unambiguous':
                                    pass
                        #print('\n'.join(sorted(map(str, items_i ^ items_j))))
                        print('--------end differing items {}-----------', which)
                        print('\n\n\n')
                    print('********* i={} j={}\n'.format(i,j))
                    show_diff('inputs')
                    show_diff('outputs')
                    print('***************************************')

#########################################################################################################################


#########################################################################################################################


# * Epilog

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    #print(_pretty_print_json(_parse_cromwell_output_str(util.file.slurp_file('/dev/shm/cromwell_out/wdl_output.txt'))))
    #print(get_docker_hash('quay.io/broadinstitute/viral-ngs'))
    #record_run_results(0,0)
    if False:
        #_copy_to_gs(git_file_path='A4.scaffolding_chosen_ref.fasta', gs_prefix='gs://sabeti-ilya-cromwell')
        _json_to_org(_dx_describe('analysis-FK7BQ580761VgPb22550gfJv'), org_file='tmp/t.org')
        #_resolve_links_in_dx_analysis(dx_analysis_id='analysis-FK7BQ580761VgPb22550gfJv', analysis_dir='exper/dxan4')
        sys.exit(0)

    if False:
        run_dx_locally(workflow_name='assemble_denovo', analysis_dxid='analysis-FJfqjg005Z3Vp5Q68jxzx5q1',
                       docker_img='quay.io/broadinstitute/viral-ngs')
    #print('\n'.join(map(str, _flatten(_json_loadf('/data/ilya-work/pipelines/an-analysis-180914-200549-b9a24c15-bb50-409e-9fa9-9996e5c5ff37-assemble_denovo/metadata.json')))))
    if False:
        compare_analysis_pairs(analysis_dirs=[os.path.join('pipelines', d) for d in os.listdir('pipelines') \
                                              if os.path.isfile(os.path.join('pipelines', d, 'metadata.json'))],
                               filter_A={'inputs.assemble_denovo.assemble.assembler': 'trinity'},
                               filter_B={'inputs.assemble_denovo.assemble.assembler': 'spades'}, metrics=()
        )
    util.cmd.main_argparse(__commands__, __doc__)

