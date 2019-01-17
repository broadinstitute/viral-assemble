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

def test_print_version():
    tools.dnanexus.DxTool().print_version()
