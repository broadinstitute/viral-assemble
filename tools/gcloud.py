'''
    Tool wrapper for accessing the Google Cloud
'''

# * imports

import logging
import collections
import os
import os.path
import subprocess
import shutil
import random
import shlex
import tempfile
import pipes
import time
import binascii
import base64

import uritools
from google.cloud import storage

import tools
import util.file
import util.misc

TOOL_NAME = 'google-cloud-storage'
TOOL_VERSION = '1.13.0'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)

# * class GCloudTool
class GCloudTool(tools.Tool):

    '''Tool wrapper for accessing the Google Cloud.
    '''

    def __init__(self, install_methods=None, use_anonymous_client=False):
        """Initialize the GCloud tool wrapper.

        Args:
          use_anonymous_client: if False (default), never use anonymous client; if True, try anonymous
            client if default auth fails; if 'always', always use anonymous client only.
        """
        if install_methods is None:
            install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION, channel='conda-forge',
                                                  executable=None)]
        tools.Tool.__init__(self, install_methods=install_methods)
        self._use_anonymous_client = 'always' if 'PYTEST_CURRENT_TEST' in os.environ else use_anonymous_client
        self._storage_client = None
        self._is_anonymous = False
        self._bucket_name_to_bucket = {}

    def version(self):
        return TOOL_VERSION
    
    def is_anonymous_client(self):
        self.get_storage_client()
        return self._is_anonymous

    def get_storage_client(self):
        if not self._storage_client:
            if self._use_anonymous_client == 'always':
                self._storage_client = storage.Client.create_anonymous_client()
                self._is_anonymous = True
            else:
                try:
                    self._storage_client = storage.Client()
                    self._is_anonymous = False
                except Exception:
                    if self._use_anonymous_client:
                        self._storage_client = storage.Client.create_anonymous_client()
                        self._is_anonymous = True
                    else:
                        raise
            util.misc.chk(self._storage_client, 'failed to initialize storage client')
            _log.info('Initialized storage client, anonymous=%s', self._is_anonymous)

        util.misc.chk(self._storage_client)
        return self._storage_client

    def get_bucket(self, bucket_name):
        """Return the bucket object for given bucket name"""
        return self.get_storage_client().bucket(bucket_name)

    def maybe_get_blob(self, gs_uri):
        """Return the blob object for given gs:// uri"""
        _log.info('getting blob from %s', gs_uri)
        gs_uri_parts = uritools.urisplit(gs_uri)
        util.misc.chk(gs_uri_parts.scheme == 'gs' and gs_uri_parts.path.startswith('/') and \
                      gs_uri_parts.query is None and gs_uri_parts.fragment is None)
        bucket = self.get_bucket(gs_uri_parts.authority)
        blob = bucket.get_blob(gs_uri_parts.path[1:])
        return blob

    def get_blob(self, gs_uri):
        return util.misc.chk(self.maybe_get_blob(gs_uri), 'could not get blob from {} {}'.format(gs_uri, gs_uri_parts.path[1:]))

    def download_object(self, gs_uri, destination_file_name):
        """Downloads a blob from the bucket."""
        blob = self.get_blob(gs_uri)
        blob.download_to_filename(destination_file_name)

    def get_metadata_for_objects(self, gs_uris, ignore_errors=False):
        """Get metadata for objects stored in Google Cloud Storage.

        TODO:
          - cache results for given gs_uri (most of our URIs don't change)
            (make sure to lock the shared cache before any update)
          - use batch request google.cloud.storage.batch.Batch
          - use 'fields' to get partial response with just the fields we want
          - parallelize with threads
        """

        def _maybe_decode(s):
            return s.decode() if hasattr(s, 'decode') else s
        
        uri2attrs = {}
        for gs_uri in gs_uris:
            blob = self.maybe_get_blob(gs_uri)
            util.misc.chk(ignore_errors or blob, 'could not get blob for {}'.format(gs_uri))
            if blob:
                md5_hex = binascii.hexlify(base64.b64decode(blob.md5_hash))
                uri2attrs[gs_uri] = collections.OrderedDict([('md5', _maybe_decode(md5_hex).upper()),
                                                             ('size', blob.size)])
        return uri2attrs

# end: class GCloudTool
