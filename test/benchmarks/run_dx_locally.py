import json
import subprocess
import os
import os.path
import shutil
import glob

import dxpy
import dxpy.bindings.dxfile_functions

import util.file

#
# given a dx analysis, run it locally
#

def get_workflow_inputs(workflow_name, t_dir):
    """Run womtool to get the inputs of the wdl workflow"""
    womtool_jar = '/data/vngs/viral-ngs-etc/mc3/envs/master_env_v25/share/cromwell/womtool.jar'
    for f in glob.glob('pipes/WDL/workflows/*.wdl'):
        shutil.copy(f, t_dir)
    for f in glob.glob('pipes/WDL/workflows/tasks/*.wdl'):
        shutil.copy(f, t_dir)
    with util.file.pushd_popd(t_dir):
        return json.loads(subprocess.check_output('java -jar ' + womtool_jar + ' inputs ' + workflow_name + '.wdl', shell=True))
    

def get_dx_val(val):
    dx_cache = '/dev/shm/dxcache'
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

def run_dx_locally(workflow_name, dxid):
    """Run a dx analysis locally"""
    analysis = dxpy.DXAnalysis(dxid=dxid)
    descr = analysis.describe()


    with util.file.tmp_dir('_workflow_copy') as t_dir:
    
        wdl_wf_inputs = get_workflow_inputs(workflow_name, t_dir)

        dx_wf_inputs = { k.split('.')[-1] : v for k, v in descr['originalInput'].items()}   # check for dups
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
            cromwell_jar = '/data/vngs/viral-ngs-etc/mc3/envs/master_env_v25/share/cromwell/cromwell.jar'
            shutil.copy(cromwell_jar, '.')
            subprocess.check_call('java -jar cromwell.jar run ' + workflow_name + '.wdl -i wf.json', shell=True)


if __name__ == '__main__':
#    print(get_workflow_inputs('assemble_denovo_with_deplete'))
    run_dx_locally(workflow_name='assemble_denovo_with_deplete', dxid='analysis-FJfqjg005Z3Vp5Q68jxzx5q1')


    

    
