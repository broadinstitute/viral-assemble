"""Unit tests for tools/docker.py"""

__author__ = "ilya@broadinstitute.org"

import logging

import pytest

import tools.docker

_log = logging.getLogger(__name__)  # pylint: disable=invalid-name

@pytest.mark.skipif(not tools.docker.DockerTool.is_docker_installed(), reason='requires docker')
def test_docker_basic():
    docker_tool = tools.docker.DockerTool()
    d_img = 'quay.io/notestaff/testing:hello-01'
    d_hash = 'sha256:cdfea997b83d847e33f1fe703cc0d37a724ad271bf5c10920831e12a6cac6528'
    assert not docker_tool.has_image_hash(d_img)
    assert docker_tool.has_image_hash(d_img+'@'+d_hash)
    assert docker_tool.get_docker_hash(d_img) == d_hash
    assert docker_tool.add_image_hash(d_img) == d_img + '@'+d_hash
    assert docker_tool.strip_image_hash(d_img + '@' + d_hash) == d_img
