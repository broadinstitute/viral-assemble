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

def test_git_annex_basic():
    ga = tools.git_annex.GitAnnexTool()
    ga.execute(['version'])

def test_git_annex_get(tmpdir_function):

    join = os.path.join
    isfile = os.path.isfile
    abspath = os.path.abspath
    relpath = os.path.relpath
    mkdir_p = util.file.mkdir_p
    pushd_popd = util.file.pushd_popd

    ga = tools.git_annex.GitAnnexTool()
    with pushd_popd(tmpdir_function):
        dir_remote = join(tmpdir_function, 'dir_remote')
        mkdir_p(dir_remote)
        mkdir_p('ga_repo')
        with pushd_popd('ga_repo'):
            ga.init_repo()
            file_A = 'testfile.txt'
            util.file.dump_file(file_A, 'some contents')
            ga.add(file_A)
            ga.commit('one file')
            assert isfile(file_A)
            assert ga._get_link_into_annex(file_A)[0] == file_A
            
            dir_remote_name = 'my_dir_remote'
            ga.initremote(dir_remote_name, 'directory', directory=dir_remote)
            
            ga.move(file_A, to_remote_name=dir_remote_name)
            assert not isfile(file_A)
            assert ga._get_link_into_annex(file_A)[0] == file_A
            ga.get(file_A)
            assert isfile(file_A)

            ga.drop(file_A)
            assert not isfile(file_A)

            file_A_abs = abspath(file_A)
            with pushd_popd('/'):
                ga.get(file_A_abs)
                assert isfile(file_A_abs)
                ga.drop(file_A_abs)
                assert not isfile(file_A_abs)
                
                file_A_rel = relpath(file_A_abs)
                ga.get(file_A_rel)
                assert isfile(file_A_abs)
                
                ga.get(file_A_rel)
                assert isfile(file_A_abs)
                ga.drop(file_A_rel)
                assert not isfile(file_A_abs)

            with pushd_popd('..'):
                file_A_link_abs = 'file_A_link_abs'
                os.symlink(file_A_abs, file_A_link_abs)
                assert not isfile(file_A_abs)
                ga.get(file_A_link_abs)
                assert isfile(file_A_abs)
                assert isfile(file_A_link_abs)
                ga.drop(file_A_link_abs)
                assert not isfile(file_A_abs)
                assert not isfile(file_A_link_abs)

                file_A_rel = relpath(file_A_abs)
                file_A_link_rel = 'file_A_link_rel'
                os.symlink(file_A_rel, file_A_link_rel)
                assert not isfile(file_A_rel)
                ga.get(file_A_link_rel)
                assert isfile(file_A_rel)
                assert isfile(file_A_link_rel)
                ga.drop(file_A_link_rel)
                assert not isfile(file_A_rel)
                assert not isfile(file_A_link_rel)



                
                
                






