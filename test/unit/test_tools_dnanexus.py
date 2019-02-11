"""Unit tests for tools/dnanexus.py"""

__author__ = "ilya@broadinstitute.org"

import os
import os.path
import sys
import collections
import argparse
import logging
import itertools

import pytest

import util.file
import util.misc
import tools.dnanexus

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

pytestmark = pytest.mark.skipif('VIRAL_NGS_DX_TESTS' not in os.environ,
                                reason='skipping DNAnexus tests because env var VIRAL_NGS_DX_TESTS is not set')

@pytest.fixture(scope='module')
def dx_tool():
    return tools.dnanexus.DxTool()

@pytest.fixture(scope='module')
def dx_file_A():
    return 'file-FVQ74G80f5zf4V5vK0GjB1q5'

def test_print_version(dx_tool):
    dx_tool.print_version()

def test_describe(dx_tool, dx_file_A):
    descr = dx_tool.describe(dx_file_A)
    assert descr['id'] == dx_file_A

def test_url_to_dxid(dx_tool, dx_file_A):
    assert dx_tool.url_to_dxid(dx_tool.DX_URI_PFX + dx_file_A) == dx_file_A




