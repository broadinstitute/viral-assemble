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

if not platform.python_version().startswith('2.7'):
    pytest.skip("skipping py27-only tests for workflow_utils", allow_module_level=True)

import util.cmd
import util.file
import util.misc
import workflow_utils
import tools.git_annex

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in workflow_utils.__commands__:
            _log.info('looking at commmand %s', cmd_name)
            parser = parser_fun(argparse.ArgumentParser())
            assert parser, 'No parser for command {}'.format(cmd_name)
            helpstring = parser.format_help()

class TestGatherFileMetadataFromAnalysisMetadata():


    def _test_fname(self, *path_elts):
        return os.path.join(util.file.get_test_input_path(), 'TestWorkflowUtils', *path_elts)

    def _call_for(self, test_data_id):
        metadata_json_fname = self._test_fname('metadata_orig.{}.json'.format(test_data_id))
        metadata_json_data = workflow_utils._json_loadf(metadata_json_fname)
        return workflow_utils._gather_file_metadata_from_analysis_metadata(metadata_json_data)

    def test_succ(self):
        """Test case of metadata from an analysis that succeeded"""
        succ_mdata = self._call_for('succ')
        succ_exp = workflow_utils._json_loadf(self._test_fname('metadata_orig.succ.exp.json'))
        assert sorted(succ_mdata.items()) == sorted(succ_exp.items())






        
        
