#!/usr/bin/env python

# * Preamble

"""
Utilities for dealing with workflows, including cloud workflows.

Commands to help jointly version code and data using tools such as git-annex and DataLad.
"""

__author__ = "ilya@broadinstitute.org"
__commands__ = []

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
import SimpleHTTPServer
import SocketServer

import dxpy
import dxpy.bindings.dxfile_functions
import boto3

import util.cmd
import util.file
import util.misc

_log = logging.getLogger(__name__)

logging.basicConfig(format="%(asctime)s - %(module)s:%(lineno)d:%(funcName)s - %(levelname)s - %(message)s")
_log.setLevel(logging.DEBUG)

# * Generic utils

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

def _is_scalar(val):
    isinstance(val, _scalar_types)

def _pretty_print_json(json_dict, **kw):
    """Return a pretty-printed version of a dict converted to json, as a string."""
    return json.dumps(json_dict, indent=4, separators=(',', ': '), **kw)

def _run(cmd):
    print('running command: ', cmd, 'cwd=', os.getcwd())
    beg_time = time.time()
    subprocess.check_call(cmd, shell=True)
    print('command succeeded in {}s: {}'.format(time.time()-beg_time, cmd))

def _run_get_output(cmd):
    print('running command: ', cmd)
    beg_time = time.time()
    output = subprocess.check_output(cmd, shell=True)
    print('command succeeded in {}s: {}'.format(time.time()-beg_time, cmd))
    return output

def _json_loads(s):
    return json.loads(s.strip(), object_pairs_hook=collections.OrderedDict)

def _json_loadf(fname):
    return _json_loads(util.file.slurp_file(fname))

def _run_get_json(cmd):
    return _json_loads(_run_get_output(cmd))

def workflow_utils_init():
    """Install the dependencies: cromwell and dxpy and git-annex."""
    _run('conda install cromwell dxpy git-annex')


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

    print('parsing val: ', val)
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
# end: def _get_dx_val(val, dx_files_dir = 'input_files')

def _parse_cromwell_output_str(cromwell_output_str):
    """Parse cromwell output"""
    assert cromwell_output_str.count('Final Outputs:') == 1
    json_beg = cromwell_output_str.index('Final Outputs:') + len('Final Outputs:')
    json_end = cromwell_output_str.index('\n}\n', json_beg) + 2
    return _json_loads(cromwell_output_str[json_beg:json_end])

# ** submit_analysis_wdl

def submit_analysis_wdl(workflow_name, analysis_inputs_from_dx_analysis, docker_img, analysis_dir, analysis_inputs_specified=None, 
                        analysis_descr='', analysis_labels=None,
                        data_repo=None, data_remote=None, cromwell_opts=''):
    """Submit a WDL analysis to a Cromwell server.

    Args:
        workflow_name: name of the workflow, from pipes/WDL/workflows
        docker_img: docker image from which all viral-ngs and WDL code is taken.
        analysis_inputs_from_dx_analysis: id of a DNAnexus analysis from which to take analysis inputs
        analysis_inputs_specified: json file specifying analysis inputs directly; overrides any given by
         analysis_inputs_from_dx_analysis
        analysis_labels: json file specifying any analysis labels
        analysis_dir: relative path where analysis results will be put.  a unique analysis ID will be appended.
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

    cromwell_tmp_dir = '/data/ilya-work/pipelines'
    util.file.mkdir_p(cromwell_tmp_dir)

    _run('git config annex.security.allowed-http-addresses "127.0.0.1 ::1 localhost"')
                # if data_remote:
                #     _run('git annex enableremote ' + data_remote)

    t_dir = os.path.join(cromwell_tmp_dir, analysis_dir + '-' + analysis_id)
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
        dx_wf_inputs = { k.split('.')[-1] : v for k, v in dx_analysis_descr['runInput'].items()}   # check for dups
        dx_wf_orig_inputs = { k.split('.')[-1] : v for k, v in dx_analysis_descr['originalInput'].items()}   # check for dups
        print('DX_WF_RUNINPUTS', '\n'.join(dx_wf_inputs.keys()))
        print('DX_WF_ORIGINPUTS', '\n'.join(dx_wf_orig_inputs.keys()))
        new_wdl_wf_inputs = collections.OrderedDict()
        for wdl_wf_input, wdl_wf_input_descr in wdl_wf_inputs.items():
            wdl_wf_input_full = wdl_wf_input
            wdl_wf_input = wdl_wf_input.split('.')[-1]

            def _set_input(dx_wf_input):
                new_wdl_wf_inputs[wdl_wf_input_full] = map(get_dx_val, dx_wf_input) \
                                                       if isinstance(dx_wf_input, list) \
                                                          else get_dx_val(dx_wf_input)

            if wdl_wf_input in analysis_inputs_specified:
                _set_input(analysis_inputs_specified[wdl_wf_input])
            elif wdl_wf_input in dx_wf_inputs:
                print('HAVE', wdl_wf_input, wdl_wf_input_descr, dx_wf_inputs[wdl_wf_input])
                _set_input(dx_wf_inputs[wdl_wf_input])
            elif '(optional' not in wdl_wf_input_descr and wdl_wf_input in dx_wf_orig_inputs:
                print('HAVE', wdl_wf_input, wdl_wf_input_descr, dx_wf_orig_inputs[wdl_wf_input])
                _set_input(dx_wf_orig_inputs[wdl_wf_input])
            else:
                print('MISSING', wdl_wf_input, wdl_wf_input_descr)
                assert '(optional' in wdl_wf_input_descr, \
                    'Missing required argument: {} {}'.format(wdl_wf_input, wdl_wf_input_descr)


        print(_pretty_print_json(new_wdl_wf_inputs))

        ################# put in the right docker ID!!  and find a place to keep the docker cache.

        with open('inputs.json', 'wt') as wf_out:
            json.dump(new_wdl_wf_inputs, wf_out, indent=4, separators=(',', ': '))

        # TODO: option to update just some of the tasks.
        # actually, when compiling WDL, should have this option -- or, actually,
        # should make a new workflow where older apps are reused for stages that have not changed.
        _run('sed -i -- "s|{}|{}|g" *.wdl'.format('quay.io/broadinstitute/viral-ngs', docker_img_hash))

        util.file.mkdir_p(output_dir)
        util.file.mkdir_p(os.path.join(output_dir, 'outputs'))
        util.file.mkdir_p(os.path.join(output_dir, 'logs'))
        util.file.mkdir_p(os.path.join(output_dir, 'call_logs'))
        #util.file.mkdir_p(os.path.join(output_dir, 'metadata'))
        wf_opts_dict = { "final_workflow_outputs_dir": os.path.join(output_dir, 'outputs'),
                         "final_workflow_log_dir": os.path.join(output_dir, 'logs'),
                         "final_call_logs_dir": os.path.join(output_dir, 'call_logs')
        }
        util.file.dump_file('cromwell_opts.json', _pretty_print_json(wf_opts_dict))
        util.file.dump_file('execution_env.json', _pretty_print_json(dict(ncpus=util.misc.available_cpu_count()),
                                                                     sort_keys=True))

        util.file.dump_file('analysis_labels.json',
                            _pretty_print_json(dict(analysis_descr=analysis_descr,
                                                    docker_img=docker_img,
                                                    docker_img_hash=docker_img_hash,
                                                    inputs_from_dx_analysis=analysis_inputs_from_dx_analysis,
                                                    analysis_id=analysis_id,
                                                    **dict(analysis_labels or {})), sort_keys=True))

        # add cromwell labels: dx project, the docker tag we ran on, etc.

        _log.info('Validating workflow')
        _run('womtool validate -i inputs.json ' + workflow_name + '.wdl')
        _log.info('Validated workflow; calling cromwell')
        _run('zip imports.zip *.wdl')
        try:
            cromwell_output_str = _run_get_output('cromwell ' + cromwell_opts + ' submit ' + workflow_name + \
                                                  '.wdl -t wdl -i inputs.json -l analysis_labels.json ' + \
                                                  ' -o cromwell_opts.json' + \
                                                  ' -p imports.zip -h http://localhost:8000')
            cromwell_returncode = 0
        except subprocess.CalledProcessError as called_process_error:
            cromwell_output_str = called_process_error.output
            cromwell_returncode = called_process_error.returncode

        _log.info('Cromwell returned with return code %d', cromwell_returncode)
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
    digest_line = _run_get_output('docker images ' + docker_img + ' --digests --no-trunc --format '
                                  '"{{.Repository}}:{{.Tag}} {{.Digest}}"')
    assert digest_line.count('\n') == 1
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
    parser.add_argument('--analysisDir', dest='analysis_dir', default='an',
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
    parser.add_argument('--cromwellOpts', dest='cromwell_opts', default='', help='command-line flags to pass to cromwell')

    util.cmd.attach_main(parser, submit_analysis_wdl, split_args=True)

__commands__.append(('submit_analysis_wdl', parser_submit_analysis_wdl))

########################################################################################################################

    

########################################################################################################################

# * RedirectToDNAnexus

class RedirectToDNAnexus(SimpleHTTPServer.SimpleHTTPRequestHandler):

   def do_GET(self):
       try:
           print('path=', self.path)
           try:
               file_info = _run_get_json('dx describe --json ' + self._get_dx_id())
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
       except (IOError, subprocess.CalledProcessError):
           self.send_error(404, 'file not found')

   def do_HEAD(self):
       self.do_GET()

   def _get_dx_id(self):
       return self.path[len('/dx/'):]

def run_dx_url_server(dummy, port=8080):
    """Start a webserver that will redirect URL requests of the form http://localhost/dx/file-xxxxxx to DNAnexus.
    This gives each dx file a stable URL, permitting such files to be added to git-annex with the addurl command."""

    server = None
    try:           
        SocketServer.TCPServer.allow_reuse_address = True
        server = SocketServer.TCPServer(("", port), RedirectToDNAnexus)
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

def parser_run_dx_url_server(parser=argparse.ArgumentParser()):
    parser.add_argument('dummy', help='Ignored argument (running command with no args just prints help)')
    parser.add_argument('--port', default=8080, help='Port on which to run the webserver')
    util.cmd.attach_main(parser, run_dx_url_server, split_args=True)

__commands__.append(('run_dx_url_server', parser_run_dx_url_server))

########################################################################################################################

# * finalize_analysis_dir

# ** CromwelServer

class CromwellServer(object):

    """Interactions with a running Cromwell server"""

    def __init__(self, host):
        self.host = host

    def _api(self, url):
        return _run_get_json('curl -X GET "http://{}/api/workflows/v1/{}" -H "accept: application/json"'.format(self.host, url))

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
            os.link(val, abspath)

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

# ** finalize_analysis_dirs impl
def finalize_analysis_dirs(cromwell_host):
    """After a submitted cromwell analysis has finished, save results to the analysis dir.
    Save metadata, mark final workflow result, make paths relative to analysis dir."""
    cromwell_server = CromwellServer(host=cromwell_host)
    for wf in cromwell_server.get_workflows():
        mdata = cromwell_server.get_metadata(wf['id'])
        assert mdata['id'] == wf['id']
        assert mdata['workflowLog'].endswith('workflow.{}.log'.format(wf['id']))
        analysis_dir = os.path.dirname(os.path.dirname(os.path.dirname(mdata['workflowLog'])))
        mdata_rel = _record_file_metadata(mdata, analysis_dir, mdata['workflowRoot'])
        mdata_fname = os.path.join(analysis_dir, 'metadata_orig.json') # mdata['workflowLog'][:-4]+'.metadata.json'
        mdata_rel_fname = os.path.join(analysis_dir, 'metadata.json') # mdata['workflowLog'][:-4]+'.metadata.json'
        if not os.path.exists(mdata_fname):
            util.file.dump_file(mdata_fname, _pretty_print_json(mdata))
            util.file.dump_file(mdata_rel_fname, _pretty_print_json(mdata_rel))
            _log.info('Wrote metadata to %s and %s', mdata_fname, mdata_rel_fname)

def parser_finalize_analysis_dirs(parser=argparse.ArgumentParser()):
    parser.add_argument('cromwell_host', default='localhost:8000', help='cromwell server hostname')
    util.cmd.attach_main(parser, finalize_analysis_dirs, split_args=True)

__commands__.append(('finalize_analysis_dirs', parser_finalize_analysis_dirs))



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
    return _json_loadf(os.path.join(analysis_dir, 'metadata.json'))

def _flatten(val, pfx=''):
    if isinstance(val, list): return functools.reduce(operator.concat,
                                                       [_flatten(v, pfx+('.' if pfx else '')+str(i)) for i, v in enumerate(val)])
    if _is_mapping(val):
        if '_is_file' in val:
            if 'md5' in val: return _flatten('md5:'+val['md5'], pfx)
            # for now, assume that if the file is from a URL then its unique id e.g. dx file-xxxx is in filename
            return _flatten(val['relpath'], pfx)
        if '_is_dir' in val:
            return _flatten(val['relpath'], pfx)
        return functools.reduce(operator.concat, [_flatten(v, pfx+('.' if pfx else '')+k) for k, v in val.items()])
    return [(pfx, val)]

def compare_analysis_pairs(analysis_dirs, filter_A, filter_B, metrics):
    """Compare pairs of analyses from `analysis_dirs` where the two analyses of a pair have
    identical values of `common_fields`, the first analysis of the pair matches the criteria in `filter_A` and
    the second matches the criteria in `filter_B`.  For such pairs, show the distribution of output 
    `metrics`.
    """

    # For each distinct combination of `common_fields` values, gather the analyses with that combination.

    flat_analyses = [ _flatten(_load_analysis_metadata(analysis_dir)) for analysis_dir in analysis_dirs if
                      os.path.isfile(os.path.join(analysis_dir, 'metadata.json')) ]
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
    for ag, a in analysis_groups.items():
        print('ag=', ag, 'a=', a)

    def _matches(a, filt):
        return all([k not in filt or filt[k] == v for k, v in a])

    found_pairs = 0
    for ag, analyses in analysis_groups.items():
        ag_A = [a for a in analyses if _matches(a, filter_A)]
        ag_B = [a for a in analyses if _matches(a, filter_B)]
        if ag_A and ag_B:
            print('\nag=', ag, '\nag_A=', ag_A[0], '\nag_B=', ag_B[0])
            found_pairs += 1
    print('found_pairs=', found_pairs)

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




# * Epilog

def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    #print(_pretty_print_json(_parse_cromwell_output_str(util.file.slurp_file('/dev/shm/cromwell_out/wdl_output.txt'))))
    #print(get_docker_hash('quay.io/broadinstitute/viral-ngs'))
    #record_run_results(0,0)
    if False:
        run_dx_locally(workflow_name='assemble_denovo', analysis_dxid='analysis-FJfqjg005Z3Vp5Q68jxzx5q1',
                       docker_img='quay.io/broadinstitute/viral-ngs')
    #print('\n'.join(map(str, _flatten(_json_loadf('/data/ilya-work/pipelines/an-analysis-180914-200549-b9a24c15-bb50-409e-9fa9-9996e5c5ff37-assemble_denovo/metadata.json')))))
    if True:
        compare_analysis_pairs(analysis_dirs=[os.path.join('pipelines', d) for d in os.listdir('pipelines') \
                                              if os.path.isfile(os.path.join('pipelines', d, 'metadata.json'))],
                               filter_A={'inputs.assemble_denovo.assemble.assembler': 'trinity'},
                               filter_B={'inputs.assemble_denovo.assemble.assembler': 'spades'}, metrics=()
        )
    #util.cmd.main_argparse(__commands__, __doc__)

