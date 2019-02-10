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

def test_print_version(dx_tool):
    dx_tool.print_version()

def test_describe(dx_tool):
    file_dxid = 'file-FVQ74G80f5zf4V5vK0GjB1q5'
    descr = dx_tool.describe(file_dxid)
    assert descr['id'] == file_dxid




