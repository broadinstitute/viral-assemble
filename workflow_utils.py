#!/usr/bin/env python
"""
Utilities for dealing with workflows, including cloud workflows.
"""

__author__ = "ilya@broadinstitute.org"
__commands__ = []

import platform
assert platform.python_version().startswith('2.7')  # dxpy requirement

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

import dxpy
import dxpy.bindings.dxfile_functions
import boto3

import util.cmd
import util.file

_log = logging.getLogger(__name__)

logging.basicConfig(format="%(asctime)s - %(module)s:%(lineno)d:%(funcName)s - %(levelname)s - %(message)s")
_log.setLevel(logging.DEBUG)

def workflow_utils_init():
    """Install the dependencies: cromwell and dxpy and git-annex."""
    subprocess.check_call('conda install cromwell dxpy git-annex')

def _pretty_print_json(json_dict):
    return json.dumps(json_dict, indent=4, separators=(',', ': '))

def _run(cmd):
    print('running command: ', cmd)
    subprocess.check_call(cmd, shell=True)
    print('command succeeded: ', cmd)

#
# given a dx analysis, run it locally
#

def get_workflow_inputs(workflow_name, t_dir, docker_img):
    """Run womtool to get the inputs of the wdl workflow"""
    with util.file.pushd_popd(t_dir):
        _run('docker run --rm ' + docker_img + ' tar cf - source/pipes/WDL > wdl.tar')
        _run('tar xvf wdl.tar')
        for f in glob.glob('source/pipes/WDL/workflows/*.wdl'):
            shutil.copy(f, '.')
        for f in glob.glob('source/pipes/WDL/workflows/tasks/*.wdl'):
            shutil.copy(f, '.')
        return json.loads(subprocess.check_output('womtool inputs ' + workflow_name + '.wdl', shell=True))
    

def _get_dx_val(val, dx_cache = '/data/dxcache'):
    """Resolve a dx value: if it is a scalar, just return that;
    if it is a dx file, fetch the file, cache it locally, and return the path to the file."""

    print('parsing val: ', val)
    util.file.mkdir_p(dx_cache)
    if isinstance(val, collections.Mapping) and '$dnanexus_link' in val:
        link = val['$dnanexus_link']
        if isinstance(link, collections.Mapping):
            return _get_dx_val(dxpy.DXAnalysis(link['analysis']).describe()['output'][link['stage']+'.'+link['field']])
        elif link.startswith('file-'):
            dxid = link
            descr = dxpy.describe(dxid)
            cache_file = os.path.join(dx_cache, dxid) + '-' + descr['name']
            file_size = int(descr['size'])
            if not os.path.isfile(cache_file) or os.path.getsize(cache_file) != file_size:
                print('fetching', dxid, 'to', cache_file)
                # TODO: check that there is enough free space (with some to spare)
                if os.path.isfile(cache_file):
                    os.unlink(cache_file)
                fs_info = os.statvfs(dx_cache)
                assert file_size < fs_info.f_bsize * fs_info.f_bavail
                dxpy.bindings.dxfile_functions.download_dxfile(dxid=dxid, filename=cache_file+'.fetching', show_progress=True)
                if os.path.getsize(cache_file+'.fetching') == file_size:
                    os.rename(cache_file+'.fetching', cache_file)
                print('fetched', dxid, 'to', cache_file)
        else:
            raise RuntimeError('Unknown dxid {}'.format(dxid))
        return cache_file
    else:
        return val

def _parse_cromwell_output_str(cromwell_output_str):
    """Parse cromwell output"""
    assert cromwell_output_str.count('Final Outputs:') == 1
    json_beg = cromwell_output_str.index('Final Outputs:') + len('Final Outputs:')
    json_end = cromwell_output_str.index('\n}\n', json_beg) + 2
    return json.loads(cromwell_output_str[json_beg:json_end])

def run_dx_locally(workflow_name, analysis_dxid, docker_img):
    """Run a dx analysis locally, using a specified version of viral-ngs.

    Notes:

       - record the actual git hash (since git tag might change)

    """
    analysis_id = create_run_id()

    _run('docker pull ' + docker_img)

    analysis = dxpy.DXAnalysis(dxid=analysis_dxid)
    analysis_descr = analysis.describe()

    with util.file.tmp_dir('_workflow_copy') as t_dir_orig:
        with util.file.pushd_popd(t_dir_orig):
            _run('git clone /ndata/cromwell-workflow-logs ' + analysis_id)
            t_dir_git = os.path.join(t_dir_orig, analysis_id)
            util.file.mkdir_p(t_dir_git)
            with util.file.pushd_popd(t_dir_git):
                _run('git annex init')
                t_dir = os.path.join(t_dir_git, analysis_id)
                util.file.mkdir_p(t_dir)
                with util.file.pushd_popd(t_dir):
                    print('TTTTTTTTTTT t_dir=', t_dir)
                    output_dir = os.path.join(t_dir, 'output')

                    wdl_wf_inputs = get_workflow_inputs(workflow_name, t_dir, docker_img=docker_img)

                    dx_wf_inputs = { k.split('.')[-1] : v for k, v in analysis_descr['runInput'].items()}   # check for dups
                    dx_wf_orig_inputs = { k.split('.')[-1] : v for k, v in analysis_descr['originalInput'].items()}   # check for dups
                    print('\n'.join(dx_wf_inputs.keys()))
                    new_wdl_wf_inputs = {}
                    for wdl_wf_input, wdl_wf_input_descr in wdl_wf_inputs.items():
                        wdl_wf_input_full = wdl_wf_input
                        wdl_wf_input = wdl_wf_input.split('.')[-1]
                        if wdl_wf_input in dx_wf_inputs:
                            print('HAVE', wdl_wf_input, wdl_wf_input_descr, dx_wf_inputs[wdl_wf_input])
                            dx_wf_input = dx_wf_inputs[wdl_wf_input]
                            new_wdl_wf_inputs[wdl_wf_input_full] = map(_get_dx_val, dx_wf_input) \
                                                                   if isinstance(dx_wf_input, list) \
                                                                      else _get_dx_val(dx_wf_input)
                        elif '(optional' not in wdl_wf_input_descr and wdl_wf_input in dx_wf_orig_inputs:
                            print('HAVE', wdl_wf_input, wdl_wf_input_descr, dx_wf_orig_inputs[wdl_wf_input])
                            dx_wf_input = dx_wf_orig_inputs[wdl_wf_input]
                            new_wdl_wf_inputs[wdl_wf_input_full] = map(_get_dx_val, dx_wf_input) \
                                                                   if isinstance(dx_wf_input, list) \
                                                                      else _get_dx_val(dx_wf_input)
                        else:
                            print('MISSING', wdl_wf_input, wdl_wf_input_descr)
                            assert '(optional' in wdl_wf_input_descr


                    print(_pretty_print_json(new_wdl_wf_inputs))

                    ################# put in the right docker ID!!  and find a place to keep the docker cache.

                    with open('wf.json', 'wt') as wf_out:
                        json.dump(new_wdl_wf_inputs, wf_out, indent=4, separators=(',', ': '))

                    if docker_img:
                        # TODO: option to update just some of the tasks.
                        # actually, when compiling WDL, should have this option -- or, actually,
                        # should make a new workflow where older apps are reused for stages that have not changed.
                        _run('sed -i -- "s|{}|{}|g" *.wdl'.format('quay.io/broadinstitute/viral-ngs', docker_img))

                    shutil.rmtree(output_dir, ignore_errors=True)
                    assert not os.path.exists(output_dir)
                    util.file.mkdir_p(output_dir)
                    util.file.mkdir_p(os.path.join(output_dir, 'outputs'))
                    util.file.mkdir_p(os.path.join(output_dir, 'logs'))
                    util.file.mkdir_p(os.path.join(output_dir, 'call_logs'))
                    #util.file.mkdir_p(os.path.join(output_dir, 'metadata'))
                    wf_opts_dict = { "final_workflow_outputs_dir": os.path.join(output_dir, 'outputs'),
                                     "final_workflow_log_dir": os.path.join(output_dir, 'logs'),
                                     "final_call_logs_dir": os.path.join(output_dir, 'call_logs')
                    }

                    with open('wf_opts.json', 'wt') as wf_opts_file:
                        json.dump(wf_opts_dict, wf_opts_file, indent=4, separators=(',', ': '))

                    # add cromwell labels: dx project, the docker tag we ran on, etc.

                    _log.info('Validating workflow')
                    _run('womtool validate -i wf.json ' + workflow_name + '.wdl')
                    _log.info('Validated workflow; calling cromwell')
                    wdl_result_str = subprocess.check_output('cromwell run ' + workflow_name + \
                                                             '.wdl -i wf.json -o wf_opts.json' + \
                                                             ' -m ' + os.path.join(output_dir, 'metadata'),
                                                             shell=True)
                    _log.info('Cromwell returned')

                    cromwell_output = _parse_cromwell_output_str(wdl_result_str)
                    #print('WDL_RESULT_STR=', wdl_result_str)
                    util.file.dump_file(os.path.join(output_dir, 'wdl_output.txt'), wdl_result_str)


                    # wdl_result = json.loads(wdl_result_str)
                    # print('---------------------------')
                    # print('WDL RESULT:')
                    # print('---------------------------')
                    # print(json.dumps(wdl_result, indent=4, separators=(',', ': ')))
                    # print('---------------------------')

                _run('sudo chown -R ec2-user * || true')
                _run('chmod -R u+w *')
                _run('git annex add')
                _run('git commit -m after_running_analysis_' + analysis_id + '.')
                _run('git annex initremote content type=directory directory=/ndata/git-annex-content/ encryption=none')
                _run('git annex copy --all --to content')
                _run('git annex dead here')
                _run('git annex sync')

def create_run_id():
    """Generate a unique ID for a run (set of steps run as part of one workflow)."""
    return util.file.string_to_file_name('-'.join(map(str, 
                                                       ('analysis', time.strftime('%Y%m%d-%H%M%S', time.localtime())[2:], 
                                                        uuid.uuid4()))))[:1024]

def get_docker_hash(docker_img):
    if docker_img.startswith('sha256:'):
        return docker_img
    digest_line = subprocess.check_output('docker images ' + docker_img + ' --digests --no-trunc --format '
                                          '"{{.Repository}}:{{.Tag}} {{.Digest}}"',
                                          shell=True)
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
def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    #print(_pretty_print_json(_parse_cromwell_output_str(util.file.slurp_file('/dev/shm/cromwell_out/wdl_output.txt'))))
    #print(get_docker_hash('quay.io/broadinstitute/viral-ngs'))
    #record_run_results(0,0)
    if True:
        run_dx_locally(workflow_name='assemble_denovo', analysis_dxid='analysis-FJfqjg005Z3Vp5Q68jxzx5q1',
                       docker_img='quay.io/broadinstitute/viral-ngs')
    #util.cmd.main_argparse(__commands__, __doc__)
