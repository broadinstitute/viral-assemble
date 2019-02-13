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

try:
    from urllib import urlencode, pathname2url
except ImportError:
    from urllib.parse import urlencode
    from urllib.request import pathname2url

import dxpy
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
    def describe(dxid):
        """Return json description for the given dxid"""
        if dxid not in DxTool._describe_result:
            DxTool._describe_result[dxid] = dxpy.describe(dxid)

        return DxTool._describe_result[dxid]
        #return _run_get_json('dx', 'describe', '--verbose', '--details', '--json', dxid)

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
