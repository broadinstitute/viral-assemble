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
join, isfile, exists, lexists, abspath, relpath = util.misc.from_module(os.path, 'join isfile exists lexists abspath relpath')

@pytest.fixture(scope='module')
def ga_tool():
    return tools.git_annex.GitAnnexTool()

def test_git_annex_basic(ga_tool):
    """Test git-annex installation"""
    ga_tool.print_version()

@pytest.fixture(scope='module')
def git_annex_repo(ga_tool, tmpdir_module):
    """Creates a git-annex repo, initializes it, makes it the current dir."""
    with pushd_popd(tmpdir_module):
        mkdir_p('ga_repo')
        with pushd_popd('ga_repo'):
            ga_tool.init_repo()
            ga_tool.initremote_external(remote_name='ldir_remote', externaltype='ldir')
            ga_tool.initremote_external(remote_name='gs_uri_remote', externaltype='gs_uri')
            ga_tool.execute_git(['config', 'annex.backend', 'MD5E'])
            ga_tool.execute_git(['config', '--type=int', 'annex.maxextensionlength', '5'])
            yield os.getcwd()

@pytest.fixture(scope='module')
def dir_remote(ga_tool, git_annex_repo, tmpdir_module):
    dir_remote_name = 'my_dir_remote_{}'.format(uuid.uuid4())
    dir_remote = join(tmpdir_module, dir_remote_name)
    mkdir_p(dir_remote)

    ga_tool.initremote(dir_remote_name, 'directory', directory=dir_remote)
    return dir_remote_name

def _a_file(sfx='', subdir=None):
    file_name = 'my_file_{}_{}.txt'.format(sfx, uuid.uuid4())
    if subdir:
        subdir = 'dir_{}_{}'.format(subdir, uuid.uuid4())
        util.file.mkdir_p(subdir)
        file_name = os.path.join(subdir, file_name)

    util.file.dump_file(file_name, file_name + 'some contents')
    return file_name

@pytest.fixture
def file_A():
    return _a_file(sfx='A')

@pytest.fixture
def file_B():
    return _a_file(sfx='B', subdir='B_dir')

def _compute_md5(fname):
    data = util.file.slurp_file(fname)
    if int(platform.python_version_tuple()[0]) > 2:
        data = data.encode('utf-8')
    return hashlib.md5(data).hexdigest().lower()

def test_git_annex_init_add_get_drop(ga_tool, git_annex_repo, dir_remote, file_A):
    """Test basic git-annex operations"""

    # TODO: add tests for when file is not in top dir of repo
    # (have fixtures for various cases)
    ga_tool.add(file_A)
    ga_tool.commit('one file')

    assert isfile(file_A)
    file_A_md5 = _compute_md5(file_A)
    file_A_size = os.path.getsize(file_A)
    file_A_key_name = file_A_md5+os.path.splitext(file_A)[1]
    file_A_key = 'MD5E-s{}--{}'.format(file_A_size, file_A_key_name)
    assert ga_tool._get_link_into_annex(file_A)[0] == file_A
    assert '.git/annex/objects/' in ga_tool._get_link_into_annex(file_A)[1]
    assert ga_tool.lookupkey(file_A) == file_A_key

    ga_tool.move(file_A, to_remote_name=dir_remote)
    assert not exists(file_A)
    assert ga_tool._get_link_into_annex(file_A)[0] == file_A
    assert '.git/annex/objects/' in ga_tool._get_link_into_annex(file_A)[1]
    assert ga_tool.lookupkey(file_A) == file_A_key

    ga_tool.get(file_A)
    assert isfile(file_A)
    assert ga_tool.is_link_into_annex(file_A)

    assert ga_tool.lookupkey(file_A) == file_A_key
    key_attrs = ga_tool.examinekey(file_A_key)
    assert key_attrs['key_name'] == file_A_key_name
    assert key_attrs['md5'] == file_A_md5
    assert key_attrs['size'] == file_A_size
    assert ga_tool.construct_key(dict(key_attrs, fname=file_A)) == file_A_key

    ga_tool.drop(file_A)
    assert not exists(file_A)
    assert ga_tool.is_link_into_annex(file_A)
    assert ga_tool.lookupkey(file_A) == file_A_key

    # test operations when the current dir is outside the repo
    file_A_abs = abspath(file_A)
    with pushd_popd('/'):
        ga_tool.get(file_A_abs)
        assert isfile(file_A_abs)
        ga_tool.drop(file_A_abs)
        assert not exists(file_A_abs)

        file_A_rel = relpath(file_A_abs)
        ga_tool.get(file_A_rel)
        assert isfile(file_A_abs)

        ga_tool.get(file_A_rel)
        assert isfile(file_A_abs)
        ga_tool.drop(file_A_rel)
        assert not exists(file_A_abs)

    with pushd_popd('..'):
        file_A_link_abs = 'file_A_link_abs'
        os.symlink(file_A_abs, file_A_link_abs)
        assert not exists(file_A_abs)
        ga_tool.get(file_A_link_abs)
        assert ga_tool.is_file_in_annex(file_A_abs)
        assert isfile(file_A_abs)
        assert isfile(file_A_link_abs)
        ga_tool.drop(file_A_link_abs)
        assert not exists(file_A_abs)
        assert not exists(file_A_link_abs)

        file_A_rel = relpath(file_A_abs)
        file_A_link_rel = 'file_A_link_rel'
        os.symlink(file_A_rel, file_A_link_rel)
        assert not exists(file_A_rel)
        ga_tool.get(file_A_link_rel)
        assert isfile(file_A_rel)
        assert isfile(file_A_link_rel)
        ga_tool.drop(file_A_link_rel)
        assert not exists(file_A_rel)
        assert not exists(file_A_link_rel)
    # end: with pushd_popd('..'):
# end: def test_git_annex_init_add_get_drop(tmpdir_function):

def test_batch_add(ga_tool, git_annex_repo, file_A, file_B):
    f2key = {}
    with ga_tool.batching() as ga_tool:
        for f in (file_A, file_B):
            assert not ga_tool.is_link_into_annex(f)
            f2key[f] = ga_tool.calckey(f)
            assert util.file.md5_for_file(f).lower() in f2key[f]
            ga_tool.add(f, now=False)
            assert not ga_tool.is_link_into_annex(f)

    for f in (file_A, file_B):
        assert ga_tool.is_link_into_annex(f)
        assert ga_tool.lookupkey(f) == f2key[f]

def test_fromkey(ga_tool, git_annex_repo, file_A, file_B):
    ga_tool.add(file_A)
    file_A_key = ga_tool.lookupkey(file_A)
    file_A_link2 = file_A+'.link2.txt'
    ga_tool.fromkey(file_A_key, file_A_link2)
    assert os.path.samefile(file_A, file_A_link2)

    with ga_tool.batching() as ga_tool:
        file_A_link3 = file_A+'.link3.txt'
        file_A_link4 = file_A+'.link4.txt'
        assert not lexists(file_A_link3)
        ga_tool.fromkey(file_A_key, file_A_link3, now=False)
        assert not lexists(file_A_link3)
        assert not lexists(file_A_link4)
        ga_tool.fromkey(file_A_key, file_A_link4, now=False)
        assert not lexists(file_A_link4)
    assert ga_tool.is_file_in_annex(file_A_link3)
    assert ga_tool.is_file_in_annex(file_A_link4)

def test_get_file_exts_for_key(ga_tool):
    test_data = {'myfile.fasta': '.fasta',
                 'myfile': '',
                 'myfile.longext': '',
                 'a': '',
                 'a.b': '.b',
                 'a.b.c': '.b.c'}
    for fname, expected in test_data.items():
        for pfx in ('', '/', '/a', '/a/b', 'a', 'b'):
            fpath = os.path.join(pfx, fname)
            assert ga_tool._get_file_exts_for_key(fpath) == expected, (fpath, expected)

def test_import_urls(ga_tool, git_annex_repo, file_A, file_B):

    gs_uris = ('gs://gcp-public-data-landsat/LC08/PRE/044/034/LC80440342016259LGN00/LC80440342016259LGN00_B1.TIF',)

    ldir_uris = (file_A, file_B)
    uris_to_import = gs_uris + ldir_uris
    #uris_to_import = ldir_uris

    url2filestat = ga_tool.import_urls(urls=uris_to_import)
    for f in uris_to_import:
        assert f in url2filestat
        assert 'git_annex_key' in url2filestat[f], 'no git_annex_key for {}: {}'.format(f, url2filestat[f])
        fn = os.path.join(str(uuid.uuid4()), os.path.basename(f))
        ga_tool.fromkey(url2filestat[f]['git_annex_key'], fn)
        ga_tool.get(fn)
        assert os.path.isfile(fn)
        assert util.file.md5_for_file(fn).lower() in url2filestat[f]['git_annex_key']
