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
import datetime
import uuid

import uritools
from google.cloud import storage
from googleapiclient import discovery
from oauth2client.client import GoogleCredentials

import tools
import util.file
import util.misc

TOOL_NAME = 'google-cloud-storage'
TOOL_VERSION = '1.13.0'

_log = logging.getLogger(__name__)
_log.setLevel(logging.DEBUG)

def _noquote(s):
    return '_noquote:' + str(s)

def _quote(s):
    s = str(s)
    if s.startswith('_noquote:'): return s[len('_noquote:'):]
    return pipes.quote(s) if hasattr(pipes, 'quote') else shlex.quote(s)

def _make_cmd(cmd, *args):
    _log.debug('ARGS=%s', args)
    return ' '.join([cmd] + [_quote(str(arg)) for arg in args if arg not in (None, '')])

def _run(cmd, *args):
    cmd = _make_cmd(cmd, *args)
    _log.info('running command: %s cwd=%s', cmd, os.getcwd())
    beg_time = time.time()
    subprocess.check_call(cmd, shell=True)
    _log.info('command succeeded in {}s: {}'.format(time.time()-beg_time, cmd))

def _run_succeeds(cmd, *args):
    try:
        _run(cmd, *args)
        return True
    except subprocess.CalledProcessError:
        return False

def _run_get_output(cmd, *args):
    cmd = _make_cmd(cmd, *args)
    _log.info('running command: %s cwd=%s', cmd, os.getcwd())
    beg_time = time.time()
    output = subprocess.check_output(cmd, shell=True)
    _log.info('command succeeded in {}s: {}'.format(time.time()-beg_time, cmd))
    return output.strip()

def _gs_stat(gs_url):
    assert gs_url.startswith('gs://')
    try:
        stat_output = _run_get_output('gsutil', 'stat', gs_url)
        result = {}
        for line in stat_output.split('\n')[1:]:
            result[line[:25].strip()] = line[25:].strip()

        if 'Hash (md5):' in result:
            result['md5'] = binascii.hexlify(result['Hash (md5):'].decode('base64'))

        return result

    except subprocess.CalledProcessError:
        return {}

def _gs_ls(gs_url):
    return _run_get_output('gsutil', 'ls', gs_url).split('\n')

def _get_gs_remote_uuid():
    return '0b42380a-45f8-4b9d-82b3-e10aaf7bab6c'  # TODO determine this from git and git-annex config

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
        return util.misc.chk(self.maybe_get_blob(gs_uri), 'could not get blob from {}'.format(gs_uri))

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

    @staticmethod
    def transfer_to_gcs(url, file_size, file_md5, bucket_name='sabeti-ilya-cromwell', project_id='viral-comp-dev'):
        """Transfer given files to GCS using Google's Data Transfer Service; return their URLs."""

        ##############

        # Create the file manifest

        with util.file.tempfname(prefix='tmp_togcs', suffix='.tsv') as manifest_fname:
            with open(manifest_fname, 'w') as manifest:
                manifest.write('TsvHttpData-1.0\n')
                manifest.write('{}\t{}\t{}\n'.format(url, file_size, file_md5))

            transfer_uuid = uuid.uuid4()
            manifest_gs_uri = 'gs://{}/transfers/{}/manifest.tsv'.format(bucket_name, transfer_uuid)
            _run('gsutil', 'cp', '-a', 'public-read', manifest_fname, manifest_gs_uri)
            assert _gs_stat(manifest_gs_uri)

        ##############

            credentials = GoogleCredentials.get_application_default()

            service = discovery.build('storagetransfer', 'v1', credentials=credentials)

            # The ID of the Google Cloud Platform Console project that the Google service
            # account is associated with.
            # Required.

            request = service.googleServiceAccounts().get(projectId=project_id)
            response = request.execute()

            # TODO: Change code below to process the `response` dict:
            _log.debug('RESPONSE TO SERVICE BUILD: %s', pformat(response))

            now = datetime.datetime.now()

            transfer_job_body = {
                'description': 'testing transfer',
                'status': 'ENABLED',
                'projectId': 'viral-comp-dev',
                'schedule': {
                    'scheduleStartDate': {
                        'day': now.day,
                        'month': now.month,
                        'year': now.year
                    },
                    'scheduleEndDate': {
                        'day': now.day,
                        'month': now.month,
                        'year': now.year
                    },
                },
                'transferSpec': {
                    'httpDataSource': {
                        'listUrl': manifest_gs_uri
                    },
                    "gcsDataSink": {
                        "bucketName": bucket_name
                    }
                }
                    # TODO: Add desired entries to the request body.
            }

            request = service.transferJobs().create(body=transfer_job_body)
            response = request.execute()

            # TODO: Change code below to process the `response` dict:
            print('RESPONSE TO TRANSFER REQUEST SUBMIT %s', pformat(response))

            url_parts = uritools.urisplit(url)
            gs_file_uri = 'gs://{}/{}'.format(bucket_name, url_parts.authority + url_parts.path)

            while True:
                _log.info('waiting for file transfer...')
                time.sleep(10)
                if _gs_stat(gs_file_uri):
                    return gs_file_uri

    # end: def transfer_to_gcs(url, file_size, file_md5, bucket_name='sabeti-ilya-cromwell', project_id='viral-comp-dev')

# end: class GCloudTool
