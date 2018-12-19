"""Unit tests for tools/gcloud.py"""

__author__ = "ilya@broadinstitute.org"

import os
import sys
import collections
import argparse
import logging
import itertools

import pytest

import util.file
import util.misc
import tools.gcloud

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

def test_gcloud_get_metadata_for_objects():
    """Test getting metadata for a set of objects"""
    gc_tool = tools.gcloud.GCloudTool(use_anonymous_client='always')
    gs_uri = 'gs://gcp-public-data-landsat/LC08/PRE/044/034/LC80440342016259LGN00/LC80440342016259LGN00_B1.TIF'
    mdata = gc_tool.get_metadata_for_objects([gs_uri])
    assert mdata[gs_uri]['size'] == 74721736
    assert mdata[gs_uri]['md5'] == 'F37E4BE81E5FAC1D33081EACDB6AF64B'




