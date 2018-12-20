"""Unit tests for tools/git_annex.py"""

__author__ = "ilya@broadinstitute.org"

import os
import sys
import collections
import argparse
import logging
import itertools
import uuid
import hashlib
import platform

import pytest

import util.file
import util.misc
import tools.git_annex

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

mkdir_p, pushd_popd = util.misc.from_module(util.file, 'mkdir_p pushd_popd')
join, isfile, abspath, relpath = util.misc.from_module(os.path, 'join isfile abspath relpath')

@pytest.fixture
def ga_tool():
    return tools.git_annex.GitAnnexTool()

@pytest.fixture(scope='module')
def ga_tool_module():
    return tools.git_annex.GitAnnexTool()

def test_git_annex_basic(ga_tool):
    """Test git-annex installation"""
    ga_tool.print_version()

@pytest.fixture(scope='module')
def git_annex_repo(ga_tool_module, tmpdir_module):
    """Creates a git-annex repo, initializes it, makes it the current dir."""
    ga_tool = ga_tool_module
    with pushd_popd(tmpdir_module):
        mkdir_p('ga_repo')
        with pushd_popd('ga_repo'):
            ga_tool.init_repo()
            ga_tool.execute_git(['config', 'annex.backend', 'MD5E'])
            ga_tool.execute_git(['config', '--type=int', 'annex.maxextensionlength', '5'])
            yield os.getcwd()

@pytest.fixture(scope='module')
def dir_remote(ga_tool_module, git_annex_repo, tmpdir_module):
    ga_tool = ga_tool_module
    dir_remote_name = 'my_dir_remote_{}'.format(uuid.uuid4())
    dir_remote = join(tmpdir_module, dir_remote_name)
    mkdir_p(dir_remote)

    ga_tool.initremote(dir_remote_name, 'directory', directory=dir_remote)
    return dir_remote_name

def _ga_file(ga_tool, sfx=''):
    file_name = 'my_file_{}_{}.txt'.format(sfx, uuid.uuid4())
    util.file.dump_file(file_name, file_name + 'some contents')
    ga_tool.add(file_name)
    ga_tool.commit('one file')
    return file_name

@pytest.fixture
def ga_file(ga_tool, git_annex_repo):
    return _ga_file(ga_tool)

@pytest.fixture(scope='module')
def file_A(ga_tool_module, git_annex_repo):
    return _ga_file(ga_tool_module, 'A')

@pytest.fixture(scope='module')
def file_B(ga_tool_module, git_annex_repo):
    return _ga_file(ga_tool_module, 'B')

def _compute_md5(fname):
    data = util.file.slurp_file(fname)
    if int(platform.python_version_tuple()[0]) > 2:
        data = data.encode('utf-8')
    return hashlib.md5(data).hexdigest().lower()

def test_git_annex_init_add_get_drop(ga_tool, git_annex_repo, dir_remote, ga_file):
    """Test basic git-annex operations"""

    # TODO: add tests for when file is not in top dir of repo
    # (have fixtures for various cases)

    assert isfile(ga_file)
    ga_file_md5 = _compute_md5(ga_file)
    ga_file_size = os.path.getsize(ga_file)
    ga_file_key_name = ga_file_md5+os.path.splitext(ga_file)[1]
    ga_file_key = 'MD5E-s{}--{}'.format(ga_file_size, ga_file_key_name)
    assert ga_tool._get_link_into_annex(ga_file)[0] == ga_file
    assert '.git/annex/objects/' in ga_tool._get_link_into_annex(ga_file)[1]
    assert ga_tool.lookupkey(ga_file) == ga_file_key

    ga_tool.move(ga_file, to_remote_name=dir_remote)
    assert not isfile(ga_file)
    assert ga_tool._get_link_into_annex(ga_file)[0] == ga_file
    assert '.git/annex/objects/' in ga_tool._get_link_into_annex(ga_file)[1]
    assert ga_tool.lookupkey(ga_file) == ga_file_key

    ga_tool.get(ga_file)
    assert isfile(ga_file)
    assert ga_tool.is_annexed_file(ga_file)

    key_attrs = ga_tool.examinekey(ga_tool.lookupkey(ga_file))
    assert key_attrs['key_name'] == ga_file_key_name
    assert key_attrs['md5'] == ga_file_md5
    assert key_attrs['size'] == ga_file_size

    ga_tool.drop(ga_file)
    assert not isfile(ga_file)
    assert ga_tool.is_annexed_file(ga_file)
    assert ga_tool.lookupkey(ga_file) == ga_file_key

    # test operations when the current dir is outside the repo
    ga_file_abs = abspath(ga_file)
    with pushd_popd('/'):
        ga_tool.get(ga_file_abs)
        assert isfile(ga_file_abs)
        ga_tool.drop(ga_file_abs)
        assert not isfile(ga_file_abs)

        ga_file_rel = relpath(ga_file_abs)
        ga_tool.get(ga_file_rel)
        assert isfile(ga_file_abs)

        ga_tool.get(ga_file_rel)
        assert isfile(ga_file_abs)
        ga_tool.drop(ga_file_rel)
        assert not isfile(ga_file_abs)

    with pushd_popd('..'):
        ga_file_link_abs = 'ga_file_link_abs'
        os.symlink(ga_file_abs, ga_file_link_abs)
        assert not isfile(ga_file_abs)
        ga_tool.get(ga_file_link_abs)
        assert isfile(ga_file_abs)
        assert isfile(ga_file_link_abs)
        ga_tool.drop(ga_file_link_abs)
        assert not isfile(ga_file_abs)
        assert not isfile(ga_file_link_abs)

        ga_file_rel = relpath(ga_file_abs)
        ga_file_link_rel = 'ga_file_link_rel'
        os.symlink(ga_file_rel, ga_file_link_rel)
        assert not isfile(ga_file_rel)
        ga_tool.get(ga_file_link_rel)
        assert isfile(ga_file_rel)
        assert isfile(ga_file_link_rel)
        ga_tool.drop(ga_file_link_rel)
        assert not isfile(ga_file_rel)
        assert not isfile(ga_file_link_rel)
    # end: with pushd_popd('..'):
# end: def test_git_annex_init_add_get_drop(tmpdir_function):

def test_batch_add(ga_tool, git_annex_repo, file_A, file_B):
    with ga_tool.batching():
        pass
    
