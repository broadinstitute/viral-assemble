"""Unit tests for workflow_utils.py"""

__author__ = "ilya@broadinstitute.org"

import os
import sys
import collections
import argparse
import logging
import itertools
import unittest

import pytest

import util.cmd
import util.file
import util.misc
import workflow_utils

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

class TestCommandHelp(unittest.TestCase):

    def test_help_parser_for_each_command(self):
        for cmd_name, parser_fun in workflow_utils.__commands__:
            _log.info('looking at commmand %s', cmd_name)
            parser = parser_fun(argparse.ArgumentParser())
            assert parser, 'No parser for command {}'.format(cmd_name)
            helpstring = parser.format_help()




