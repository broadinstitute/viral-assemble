"""Unit tests for tools/gcloud.py"""

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
import tools.gcloud

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

def test_gcloud_get_metadata_for_objects():
    """Test getting metadata for a set of objects"""
    gcloud_tool = tools.gcloud.GCloudTool(use_anonymous_client='always')
    gs_uri = 'gs://gcp-public-data-landsat/LC08/PRE/044/034/LC80440342016259LGN00/LC80440342016259LGN00_B1.TIF'
    mdata = gcloud_tool.get_metadata_for_objects([gs_uri])
    assert mdata[gs_uri]['size'] == 74721736
    assert mdata[gs_uri]['md5'] == 'F37E4BE81E5FAC1D33081EACDB6AF64B'

def test_download_object(tmpdir_function):
    gcloud_tool = tools.gcloud.GCloudTool(use_anonymous_client='always')
    gs_uri = 'gs://gcp-public-data-landsat/LC08/PRE/044/034/LC80440342016259LGN00/LC80440342016259LGN00_B1.TIF'
    tmp_fname = os.path.join(tmpdir_function, 'tempfile.tif')
    gcloud_tool.download_object(gs_uri=gs_uri, destination_file_name=tmp_fname)
    assert os.path.isfile(tmp_fname)
    assert os.path.getsize(tmp_fname) == 74721736
    assert util.file.md5_for_file(tmp_fname) == 'F37E4BE81E5FAC1D33081EACDB6AF64B'

