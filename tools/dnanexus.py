'''
    Tool wrapper for the DNAnexus CLI.
'''

# * imports

import logging
import os
import os.path
import subprocess
import shutil
import random
import shlex
import tempfile
import pipes
import time
import platform
import functools
import collections
import copy
import json

try:
    from urllib import urlencode, pathname2url
except ImportError:
    from urllib.parse import urlencode
    from urllib.request import pathname2url

import dxpy
import dxpy.bindings.dxproject
import dxpy.exceptions
import uritools

import tools
import util.file
import util.misc

TOOL_NAME = 'dxpy'
TOOL_VERSION = '0.273.0'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)

# * class DxTool

class DxTool(tools.Tool):

    '''Tool wrapper for DNAnexus CLI'''

# ** init, execute
    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION,
                                                  executable='dx')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, args):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path()] + list(map(str, args))
        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    def print_version(self):
        """Print dx cli version"""
        self.execute(['--version'])

    DX_URI_PFX='dx://'

    @staticmethod
    def url_to_dxid(url):
        """Given a dx://file-xxx URL, return the dxid in it"""
        return uritools.urisplit(url).authority

    @staticmethod
    def dxid_to_url(dxid, include_filename=False):
        """Given a dxid, return a dx://file-xxx URL for it"""
        url = DxTool.DX_URI_PFX + dxid
        if include_filename:
            url += '/' + DxTool.describe(dxid)['name']
        return url
    
    # @classmethod
    # def dx_make_download_url(cls, dxid, duration='2h'):
    #     return _run_get_output('dx', 'make_download_url', dxid, '--duration', duration)

    # @classmethod
    # def standardize_dx_url(cls, url):
    #     dxid = _url_to_dxid(url)
    #     dx_descr = _dx_describe(_url_to_dxid(url))
    #     return cls.DX_URI_PFX + dxid + '/' + pathname2url(dx_descr['name'])

#    @util.misc.memoize_persist(to_picklable=functools.partial(json.dumps, separators=(',',':')),
#                               from_picklable=_json_loads)

    _describe_result = {}


    @staticmethod
    def describe(dxid_or_link):
        """Return json description for the given dxid"""
        key = json.dumps(dxid_or_link, sort_keys=True)
        if key not in DxTool._describe_result:
            DxTool._describe_result[key] = dxpy.describe(dxid_or_link)

        return DxTool._describe_result[key]
        #return _run_get_json('dx', 'describe', '--verbose', '--details', '--json', dxid)

    VIRAL_NGS_DX_CI_PROJECT = 'project-F8PQ6380xf5bK0Qk0YPjB17P' # viral-ngs public CI project
    

    def determine_viral_ngs_dx_analysis_docker_img(self, mdata):
        """Determine the viral-ngs docker image used for a given dx analysis"""
        # determine the viral-ngs version
        executable_descr = self.describe(dxpy.dxlink(object_id=mdata['stages'][0]['execution']['executable'],
                                                            project_id=self.VIRAL_NGS_DX_CI_PROJECT))
        util.misc.chk('dxWDL' in executable_descr['tags'])
        util.misc.chk(executable_descr['createdBy'] == {'user': 'user-sabeti_ci'})
        util.misc.chk(executable_descr['folder'].startswith('/build/quay.io/broadinstitute/viral-ngs'))

        exe_folder_parts = executable_descr['folder'].split('/')
        docker_img = '/'.join(exe_folder_parts[2:5]) + ':' + exe_folder_parts[5]

        # determine an update to runInputs
        ci_proj = dxpy.dxproject.DXProject(self.VIRAL_NGS_DX_CI_PROJECT)
        orig_ci_workflow = [item['id'] for item in ci_proj.list_folder(executable_descr['folder'])['objects']
                            if item['id'].startswith('workflow')][0]
        orig_workflow_descr = self.describe(orig_ci_workflow)
        mdata['workflowName'] = orig_workflow_descr['name']
        this_workflow_descr = mdata['workflow']
        util.misc.chk(len(this_workflow_descr['stages']) == len(orig_workflow_descr['stages']))
        for this_stage, orig_stage in zip(this_workflow_descr['stages'], orig_workflow_descr['stages']):
            util.misc.chk(this_stage['id'] == orig_stage['id'])
            util.misc.chk(this_stage['name'] == orig_stage['name'])
            util.misc.chk(this_stage['executable'] == orig_stage['executable'])
            
        util.misc.chk(len(this_workflow_descr['inputSpec']) == len(orig_workflow_descr['inputSpec']))
        for this_input_spec, orig_input_spec in zip(this_workflow_descr['inputSpec'], orig_workflow_descr['inputSpec']):
            util.misc.chk(this_input_spec['name'] == orig_input_spec['name'])
            util.misc.chk(this_input_spec['class'] == orig_input_spec['class'])
            util.misc.chk(this_input_spec['group'] == orig_input_spec['group'])
            if 'default' in this_input_spec and this_input_spec['default'] != orig_input_spec.get('default', None) and \
               this_input_spec['name'] in mdata['input'] and \
               this_input_spec['name'] not in mdata['runInput']:
                mdata['runInput'][this_input_spec['name']] = copy.deepcopy(mdata['input'][this_input_spec['name']])

        return docker_img

    @staticmethod
    def resolve_dx_link_to_dx_file_id_or_value(val, dx_analysis_id):
        """Resolve DNAnexus links, including indirect ones that point to the output of an analysis stage,
        to either direct link to a DNAnexus file or a simple value.
        If `val` represents a DNAnexus link to a file, return {$dnanexus_link: file-xxxx}.
        If `val` represents a DNAnexus link to a value, return that value.
        Else, return `val` unchanged.
        """

        _maps = util.misc.maps
        _is_str = util.misc.is_str

        #print('parsing val: ', val)
        recurse = functools.partial(DxTool.resolve_dx_link_to_dx_file_id_or_value,
                                    dx_analysis_id=dx_analysis_id)
        if not _maps(val, '$dnanexus_link'):
            return val

        link = val['$dnanexus_link']
        if _maps(link, 'stage') and ('field' in link or 'outputField' in link):
            _log.debug('link is %s', link)
            linked_analysis_descr = DxTool.describe(link.get('analysis', dx_analysis_id))
            linked_field = link['field'] if 'field' in link else link['outputField']
            if link['stage']+'.'+linked_field not in linked_analysis_descr['output']:
                # this can happen if the stage has failed.
                return {'$dnanexus_failed_stage_output': None}
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

    @staticmethod
    def resolve_dx_links_in_dx_analysis_descr(dx_analysis_descr):
        """Resolve dx links in a dx analysis description"""
        print('TYPE IS', type(dx_analysis_descr))
        resolved = util.misc.transform_json_data(dx_analysis_descr,
                                                 functools.partial(DxTool.resolve_dx_link_to_dx_file_id_or_value,
                                                                   dx_analysis_id=dx_analysis_descr['id']))
        def _dx_link_to_url(val):
            if dxpy.is_dxlink(val):
                util.misc.chk(val['$dnanexus_link'].startswith('file-'))
                return DxTool.dxid_to_url(val['$dnanexus_link'], include_filename=True)
            else:
                return val
        resolved = util.misc.transform_json_data(resolved, _dx_link_to_url)
        return resolved

# end: class DxTool(tools.Tool)
