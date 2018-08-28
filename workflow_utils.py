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

import dxpy
import dxpy.bindings.dxfile_functions

import util.cmd
import util.file

log = logging.getLogger(__name__)

def workflow_utils_init():
    """Install the dependencies"""
    subprocess.check_call('conda install cromwell dxpy')

#
# given a dx analysis, run it locally
#

def get_workflow_inputs(workflow_name, t_dir):
    """Run womtool to get the inputs of the wdl workflow"""
    for f in glob.glob('pipes/WDL/workflows/*.wdl'):
        shutil.copy(f, t_dir)
    for f in glob.glob('pipes/WDL/workflows/tasks/*.wdl'):
        shutil.copy(f, t_dir)
    with util.file.pushd_popd(t_dir):
        return json.loads(subprocess.check_output('womtool inputs ' + workflow_name + '.wdl', shell=True))
    

def get_dx_val(val):
    dx_cache = '/dev/shm/dxcache'
    util.file.mkdir_p(dx_cache)
    if '$dnanexus_link' in val:
        dxid = val['$dnanexus_link']
        assert dxid.startswith('file-')
        descr = dxpy.describe(dxid)
        cache_file = os.path.join(dx_cache, dxid) + '-' + descr['name']
        file_size = int(descr['size'])
        if not os.path.isfile(cache_file) or os.path.getsize(cache_file) != file_size:
            print('fetching', dxid, 'to', cache_file)
            if os.path.isfile(cache_file):
                os.unlink(cache_file)
            dxpy.bindings.dxfile_functions.download_dxfile(dxid=dxid, filename=cache_file+'.fetching', show_progress=True)
            if os.path.getsize(cache_file+'.fetching') == file_size:
                os.rename(cache_file+'.fetching', cache_file)
            print('fetched', dxid, 'to', cache_file)
        return cache_file
    else:
        return val

def run_dx_locally(workflow_name, analysis_dxid, output_dir, docker_tag=None):
    """Run a dx analysis locally"""
    analysis = dxpy.DXAnalysis(dxid=analysis_dxid)
    analysis_descr = analysis.describe()

    with util.file.tmp_dir('_workflow_copy') as t_dir:
    
        wdl_wf_inputs = get_workflow_inputs(workflow_name, t_dir)

        dx_wf_inputs = { k.split('.')[-1] : v for k, v in analysis_descr['originalInput'].items()}   # check for dups
        print('\n'.join(dx_wf_inputs.keys()))
        new_wdl_wf_inputs = {}
        for wdl_wf_input, wdl_wf_input_descr in wdl_wf_inputs.items():
            wdl_wf_input_full = wdl_wf_input
            wdl_wf_input = wdl_wf_input.split('.')[-1]
            if wdl_wf_input in dx_wf_inputs:
                print('HAVE', wdl_wf_input, wdl_wf_input_descr, dx_wf_inputs[wdl_wf_input])
                dx_wf_input = dx_wf_inputs[wdl_wf_input]
                new_wdl_wf_inputs[wdl_wf_input_full] = map(get_dx_val, dx_wf_input) if isinstance(dx_wf_input, list) else get_dx_val(dx_wf_input)
            else:
                #print('MISSING', wdl_wf_input, wdl_wf_input_descr)
                assert '(optional' in wdl_wf_input_descr


        print(json.dumps(new_wdl_wf_inputs, indent=4, separators=(',', ': ')))

        ################# put in the right docker ID!!  and find a place to keep the docker cache.

        with util.file.pushd_popd(t_dir):
            with open('wf.json', 'wt') as wf_out:
                json.dump(new_wdl_wf_inputs, wf_out, indent=4, separators=(',', ': '))

            if docker_tag:
                # TODO: option to update just some of the tasks.
                # actually, when compiling WDL, should have this option -- or, actually,
                # should make a new workflow where older apps are reused for stages that have not changed.
                subprocess.check_call('sed -i -- "s|{}|{}|g" *.wdl'.format('quay.io/broadinstitute/viral-ngs', docker_tag),
                                      shell=True)

            util.file.mkdir_p(output_dir)
            util.file.mkdir_p(os.path.join(output_dir, 'outputs'))
            util.file.mkdir_p(os.path.join(output_dir, 'logs'))
            util.file.mkdir_p(os.path.join(output_dir, 'call_logs'))
            wf_opts_dict = { "final_workflow_outputs_dir": os.path.join(output_dir, 'outputs'),
                             "final_workflow_log_dir": os.path.join(output_dir, 'logs'),
                             "final_call_logs_dir": os.path.join(output_dir, 'call_logs')
            }
                             
            with open('wf_opts.json', 'wt') as wf_opts_file:
                json.dump(wf_opts_dict, wf_opts_file, indent=4, separators=(',', ': '))

            wdl_result_str = subprocess.check_output('cromwell run ' + workflow_name + '.wdl -i wf.json -o wf_opts.json', shell=True)
            print('WDL_RESULT_STR=', wdl_result_str)
            util.file.dump_file(os.path.join(output_dir, 'wdl_output.txt'), wdl_result_str)
            
                                      
            # wdl_result = json.loads(wdl_result_str)
            # print('---------------------------')
            # print('WDL RESULT:')
            # print('---------------------------')
            # print(json.dumps(wdl_result, indent=4, separators=(',', ': ')))
            # print('---------------------------')


# =======================
def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    run_dx_locally(workflow_name='assemble_denovo_with_deplete', analysis_dxid='analysis-FJfqjg005Z3Vp5Q68jxzx5q1',
                   docker_tag='quay.io/broadinstitute/viral-ngs-dev:1.21.2-5-g394c91e-is-1808241104-refine-assembly-add-reporting',
                   output_dir='/dev/shm/cromwell_out')
    #util.cmd.main_argparse(__commands__, __doc__)
