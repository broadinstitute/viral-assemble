import dxpy
import json
import subprocess
import util.file
import shutil
import glob

#
# given a dx analysis, run it locally
#

def get_workflow_inputs(workflow_name):
    """Run womtool to get the inputs of the wdl workflow"""
    womtool_jar = '/data/vngs/viral-ngs-etc/mc3/envs/master_env_v25/share/cromwell/womtool.jar'
    with util.file.tmp_dir('_get_workflow_inputs') as t_dir:
        for f in glob.glob('pipes/WDL/workflows/*.wdl'):
            shutil.copy(f, t_dir)
        for f in glob.glob('pipes/WDL/workflows/tasks/*.wdl'):
            shutil.copy(f, t_dir)
        with util.file.pushd_popd(t_dir):
            return json.loads(subprocess.check_output('java -jar ' + womtool_jar + ' inputs ' + workflow_name + '.wdl', shell=True))
    

def run_dx_locally(workflow_name, dxid):
    """Run a dx analysis locally"""
    analysis = dxpy.DXAnalysis(dxid=dxid)
    descr = analysis.describe()

    wdl_wf_inputs = get_workflow_inputs(workflow_name)
    
    dx_wf_inputs = { k.split('.')[-1] : v for k, v in descr['originalInput'].items()}   # check for dups
    print('\n'.join(dx_wf_inputs.keys()))
    new_wdl_wf_inputs = {}
    for wdl_wf_input, wdl_wf_input_descr in wdl_wf_inputs.items():
        wdl_wf_input = wdl_wf_input.split('.')[-1]
        if wdl_wf_input in dx_wf_inputs:
            print('HAVE', wdl_wf_input, wdl_wf_input_descr, dx_wf_inputs[wdl_wf_input])
            
        else:
            #print('MISSING', wdl_wf_input, wdl_wf_input_descr)
            assert '(optional' in wdl_wf_input_descr

if __name__ == '__main__':
#    print(get_workflow_inputs('assemble_denovo_with_deplete'))
    run_dx_locally(workflow_name='assemble_denovo_with_deplete', dxid='analysis-FJfqjg005Z3Vp5Q68jxzx5q1')


    

    
