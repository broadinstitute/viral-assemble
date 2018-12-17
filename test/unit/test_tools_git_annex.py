"""Unit tests for tools/git_annex.py"""

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
import tools.git_annex

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

def test_git_annex_basic():
    """Test git-annex installation"""
    ga = tools.git_annex.GitAnnexTool()
    ga.print_version()

def test_git_annex_init_add_get_drop(tmpdir_function):
    """Test basic git-annex operations"""

    join, isfile, abspath, relpath = util.misc.from_module(os.path, 'join isfile abspath relpath')
    mkdir_p, pushd_popd = util.misc.from_module(util.file, 'mkdir_p pushd_popd')

    ga = tools.git_annex.GitAnnexTool()
    with pushd_popd(tmpdir_function):
        dir_remote = join(tmpdir_function, 'dir_remote')
        mkdir_p(dir_remote)
        mkdir_p('ga_repo')
        with pushd_popd('ga_repo'):
            ga.init_repo()
            ga.execute_git(['config', 'annex.backend', 'MD5E'])
            ga.execute_git(['config', '--type=int', 'annex.maxextensionlength', '5'])

            file_A = 'testfile.txt'
            util.file.dump_file(file_A, 'some contents')
            ga.add(file_A)
            ga.commit('one file')
            assert isfile(file_A)
            assert ga._get_link_into_annex(file_A)[0] == file_A
            assert '.git/annex/objects/' in ga._get_link_into_annex(file_A)[1]
            assert ga.lookupkey(file_A) == 'MD5E-s13--220c7810f41695d9a87d70b68ccf2aeb.txt'
            
            dir_remote_name = 'my_dir_remote'
            ga.initremote(dir_remote_name, 'directory', directory=dir_remote)
            
            ga.move(file_A, to_remote_name=dir_remote_name)
            assert not isfile(file_A)
            assert ga._get_link_into_annex(file_A)[0] == file_A
            assert '.git/annex/objects/' in ga._get_link_into_annex(file_A)[1]
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
            # end: with pushd_popd('..'):
        # end: with pushd_popd('ga_repo'):
    # end: with pushd_popd(tmpdir_function):
# end: def test_git_annex_init_add_get_drop(tmpdir_function):
