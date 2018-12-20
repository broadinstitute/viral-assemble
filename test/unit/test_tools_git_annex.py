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
ga = tools.git_annex.GitAnnexTool()

def test_git_annex_basic():
    """Test git-annex installation"""
    ga.print_version()

@pytest.fixture(scope='module')
def git_annex_repo(tmpdir_module):
    """Creates a git-annex repo, initializes it, makes it the current dir."""
    with pushd_popd(tmpdir_module):
        mkdir_p('ga_repo')
        with pushd_popd('ga_repo'):
            ga.init_repo()
            ga.execute_git(['config', 'annex.backend', 'MD5E'])
            ga.execute_git(['config', '--type=int', 'annex.maxextensionlength', '5'])
            yield os.getcwd()

@pytest.fixture(scope='module')
def dir_remote(git_annex_repo, tmpdir_module):
    dir_remote_name = 'my_dir_remote_{}'.format(uuid.uuid4())
    dir_remote = join(tmpdir_module, dir_remote_name)
    mkdir_p(dir_remote)

    ga.initremote(dir_remote_name, 'directory', directory=dir_remote)
    yield dir_remote_name

def _ga_file(sfx=''):
    file_name = 'my_file_{}_{}.txt'.format(sfx, uuid.uuid4())
    util.file.dump_file(file_name, file_name + 'some contents')
    ga.add(file_name)
    ga.commit('one file')
    return file_name

@pytest.fixture
def ga_file(git_annex_repo):
    yield _ga_file()

@pytest.fixture(scope='module')
def file_A(git_annex_repo):
    yield _ga_file('A')

@pytest.fixture(scope='module')
def file_B(git_annex_repo):
    yield _ga_file('B')

def _compute_md5(fname):
    data = util.file.slurp_file(fname)
    if platform.python_version_tuple()[0] > 2:
        data = data.encode('utf-8')
    return hashlib.md5(data).hexdigest().lower()

def test_git_annex_init_add_get_drop(git_annex_repo, dir_remote, ga_file):
    """Test basic git-annex operations"""

    assert isfile(ga_file)
    ga_file_md5 = _compute_md5(ga_file)
    ga_file_size = os.path.getsize(ga_file)
    ga_file_key_name = ga_file_md5+os.path.splitext(ga_file)[1]
    ga_file_key = 'MD5E-s{}--{}'.format(ga_file_size, ga_file_key_name)
    assert ga._get_link_into_annex(ga_file)[0] == ga_file
    assert '.git/annex/objects/' in ga._get_link_into_annex(ga_file)[1]
    assert ga.lookupkey(ga_file) == ga_file_key

    ga.move(ga_file, to_remote_name=dir_remote)
    assert not isfile(ga_file)
    assert ga._get_link_into_annex(ga_file)[0] == ga_file
    assert '.git/annex/objects/' in ga._get_link_into_annex(ga_file)[1]
    assert ga.lookupkey(ga_file) == ga_file_key

    ga.get(ga_file)
    assert isfile(ga_file)
    assert ga.is_annexed_file(ga_file)

    key_attrs = ga.examinekey(ga.lookupkey(ga_file))
    assert key_attrs['key_name'] == ga_file_key_name
    assert key_attrs['md5'] == ga_file_md5
    assert key_attrs['size'] == ga_file_size

    ga.drop(ga_file)
    assert not isfile(ga_file)
    assert ga.is_annexed_file(ga_file)
    assert ga.lookupkey(ga_file) == ga_file_key

    # test operations when the current dir is outside the repo
    ga_file_abs = abspath(ga_file)
    with pushd_popd('/'):
        ga.get(ga_file_abs)
        assert isfile(ga_file_abs)
        ga.drop(ga_file_abs)
        assert not isfile(ga_file_abs)

        ga_file_rel = relpath(ga_file_abs)
        ga.get(ga_file_rel)
        assert isfile(ga_file_abs)

        ga.get(ga_file_rel)
        assert isfile(ga_file_abs)
        ga.drop(ga_file_rel)
        assert not isfile(ga_file_abs)

    with pushd_popd('..'):
        ga_file_link_abs = 'ga_file_link_abs'
        os.symlink(ga_file_abs, ga_file_link_abs)
        assert not isfile(ga_file_abs)
        ga.get(ga_file_link_abs)
        assert isfile(ga_file_abs)
        assert isfile(ga_file_link_abs)
        ga.drop(ga_file_link_abs)
        assert not isfile(ga_file_abs)
        assert not isfile(ga_file_link_abs)

        ga_file_rel = relpath(ga_file_abs)
        ga_file_link_rel = 'ga_file_link_rel'
        os.symlink(ga_file_rel, ga_file_link_rel)
        assert not isfile(ga_file_rel)
        ga.get(ga_file_link_rel)
        assert isfile(ga_file_rel)
        assert isfile(ga_file_link_rel)
        ga.drop(ga_file_link_rel)
        assert not isfile(ga_file_rel)
        assert not isfile(ga_file_link_rel)
    # end: with pushd_popd('..'):
# end: def test_git_annex_init_add_get_drop(tmpdir_function):
