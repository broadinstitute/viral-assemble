"""Unit tests for workflow_utils.py"""

__author__ = "ilya@broadinstitute.org"

import os
import sys
import collections
import argparse
import logging
import itertools
import platform
import unittest

import pytest

import util.cmd
import util.file
import util.misc
import workflow_utils
import tools.git_annex
import tools.cromwell

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in workflow_utils.__commands__:
            _log.info('looking at commmand %s', cmd_name)
            parser = parser_fun(argparse.ArgumentParser())
            assert parser, 'No parser for command {}'.format(cmd_name)
            helpstring = parser.format_help()


@pytest.fixture(scope='module')
def gcloud_tool():
    return tools.gcloud.GCloudTool(use_anonymous_client='always')

class TestGatherFileMetadataFromAnalysisMetadata():

    def _test_fname(self, *path_elts):
        return os.path.join(util.file.get_test_input_path(), 'TestWorkflowUtils', *path_elts)

    def _call_for(self, test_data_id, gcloud_tool):
        metadata_json_fname = self._test_fname('metadata_orig.{}.json'.format(test_data_id))
        metadata_json_data = workflow_utils._json_loadf(metadata_json_fname)
        assert gcloud_tool.is_anonymous_client()
        return workflow_utils._gather_file_metadata_from_analysis_metadata(metadata_json_data, gcloud_tool=gcloud_tool)

    def test_succ(self, gcloud_tool):
        """Test case of metadata from an analysis that succeeded"""
        succ_mdata = self._call_for('succ', gcloud_tool)
        #workflow_utils._write_json('/tmp/cmp.json', **succ_mdata)
        succ_exp = workflow_utils._json_loadf(self._test_fname('metadata_orig.succ.exp.json'))
        for file_name, file_mdata in succ_exp.items():
            assert set(file_mdata.items()) <= set(succ_mdata[file_name].items())

    def test_callCaching_crc32(self, gcloud_tool):
        """Test case of metadata from an analysis where some callCaching info is recorded with crc32 rather than md5"""
        callCaching_crc32_mdata = self._call_for('callCaching_crc32', gcloud_tool)
        callCaching_crc32_mdata_exp = workflow_utils._json_loadf(self._test_fname('metadata_orig.callCaching_crc32.exp.json'))

        for file_name, file_mdata in callCaching_crc32_mdata_exp.items():
            assert set(file_mdata.items()) <= set(callCaching_crc32_mdata[file_name].items())

# end: class TestGatherFileMetadataFromAnalysisMetadata()

@pytest.fixture(scope='module')
def cromwell_server():
    """Runs a Cromwell server"""
    with tools.cromwell.CromwellTool().cromwell_server() as server:
        yield server

def test_starting_cromwell_server(cromwell_server):
    _log.info('SERVER HEALTH IS %s', cromwell_server.health())

