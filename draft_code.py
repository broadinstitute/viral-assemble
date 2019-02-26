# * Various bits of code that we may want to reuse, but do not use right now

# ** from workflow_utils

#                 _run('sudo chown -R $USER . || true')
#                 _run('git annex add')
#                 _run('git commit -m after_running_analysis_' + analysis_id + '.')
# #                _run('git annex initremote content type=directory directory=/ndata/git-annex-content/ encryption=none')
#                 _run('git annex move --all --to {} -J{}'.format(data_remote or 'origin', util.misc.available_cpu_count()))
#                 _run('git annex dead here')
#                 # pull and merge here first? and try the sync several times, after a pause maybe
#                 _run('git annex sync --message git_annex_sync_analysis_{}'.format(analysis_id))

#                 # enable cleanup
#                 _run('chmod -R u+w . || true')


#        import sys
#        sys.exit(0)


        # util.file.dump_file(os.path.join(output_dir, 'cromwell_output.txt'), cromwell_output_str)
        # _run('sed -i -- "s|{}|{}|g" cromwell_execution_metadata.json'.format(t_dir+'/', ''))

        # if cromwell_returncode == 0:
        #     def make_paths_relative(v):
        #         if _is_str(v) and os.path.isabs(v) and v.startswith(t_dir):
        #             return os.path.relpath(v, t_dir)
        #         if isinstance(v, list):
        #             return list(map(make_paths_relative, v))
        #         return v
        #     cromwell_output_json = {k: make_paths_relative(v)
        #                             for k, v in _parse_cromwell_output_str(cromwell_output_str).items()}
        #     util.file.dump_file('outputs.json',
        #                         _pretty_print_json(cromwell_output_json))
        #     util.file.make_empty('analysis_succeeded.txt')
        # else:
        #     util.file.make_empty('analysis_failed.txt')

        # def _rewrite_paths_to_gs(val):
        #     if _is_mapping(val): return {k: _rewrite_paths_to_gs(v) for k, v in val.items()}
        #     if isinstance(val, list): return list(map(_rewrite_paths_to_gs, val))
        #     if _is_str(val) and os.path.isfile(val):
        #         assert val.startswith(input_files_dir)
        #         val = 'gs://sabeti-ilya-cromwell/cromwell-inputs/' + analysis_id + '/' + val[len(input_files_dir)+1:]
        #     return val
        # #_run('gsutil cp -r ' + input_files_dir + '/* gs://sabeti-ilya-cromwell/cromwell-inputs/' + analysis_id + '/')

        # _write_json('inputs.json', **new_wdl_wf_inputs)
#            json.dump(_rewrite_paths_to_gs(new_wdl_wf_inputs), wf_out, indent=4, separators=(',', ': '))

def _resolve_file_values(input_name, val, input_files_dir):
    """For workflow inputs that point to a file, ensure that the file is present in the filesystem in which the 
    analysis will be run.

    Args:
      input_name: the name of an input to a workflow
      val: the value of the input.  May be a scalar, or a composite value such as a list or a map.
           May be a value representing a file.
    """

    #Resolve a dx value: if it is a scalar, just return that;
    #if it is a dx file, fetch the file, cache it locally, and return the path to the file.

    #print('parsing val: ', val)
    util.file.mkdir_p(dx_files_dir)
    if isinstance(val, collections.Mapping) and '$dnanexus_link' in val:
        link = val['$dnanexus_link']
        if isinstance(link, collections.Mapping) and 'analysis' in link:
            _log.debug('link is %s', link)
            return _get_dx_val(dxpy.DXAnalysis(link['analysis']).describe()['output'][link['stage']+'.'+link['field']],
                               dx_files_dir)
        elif (_is_str(link) and link.startswith('file-')) or (isinstance(link, collections.Mapping) and 'id' in link and _is_str(link['id']) and link['id'].startswith('file-')):
            if _is_str(link) and link.startswith('file-'):
                dxid = link
            else:
                dxid = link['id']
            descr = dxpy.describe(dxid)
            dx_file = os.path.join(dx_files_dir, dxid) + '-' + descr['name']
            file_size = int(descr['size'])

            # see if the file is cached in git-annex
            ga_mdata = _run_get_json('git annex metadata --json --key=WORM-s0-m0--dx-' + dxid)
            if ga_mdata['fields']:
                ga_key = ga_mdata['fields']['ga_key'][0]
                _run('git annex get --key ' + ga_key)
                _run('git annex fromkey ' + ga_key + ' ' + dx_file)

            if not os.path.isfile(dx_file) or os.path.getsize(dx_file) != file_size:
                _log.debug('fetching %s to %s', dxid, dx_file)
                # TODO: check that there is enough free space (with some to spare)
                if os.path.isfile(dx_file):
                    os.unlink(dx_file)
                fs_info = os.statvfs(dx_files_dir)
                assert file_size < fs_info.f_bsize * fs_info.f_bavail
                dxpy.bindings.dxfile_functions.download_dxfile(dxid=dxid, filename=dx_file+'.fetching', show_progress=True)
                assert os.path.getsize(dx_file+'.fetching') == file_size
                os.rename(dx_file+'.fetching', dx_file)
                _log.debug('fetched %s to %s', dxid, dx_file)
                _log.debug('curdir is %s', os.getcwd())
                _run('git annex add ' + dx_file)
                ga_key = _run_get_output('git annex lookupkey ' + dx_file).strip()
                _run('git annex metadata -s dxid+=' + dxid + ' ' + dx_file)
                # record a mapping from the dxid to the git-annex key
                _run('git annex metadata --key=WORM-s0-m0--{} -s ga_key={}'.format('dx-'+dxid, ga_key))
                # register a URL that can be used to re-fetch this file from DNAnexus;
                # the URL is only valid if the 'run_dx_url_server' command is running.
                _run('git annex registerurl ' + ga_key + ' ' + ' http://localhost:8080/dx/' + dxid)
        else:
            raise RuntimeError('Cannot parse dx link {}'.format(link))
        return dx_file
    else:
        return val
# end: def _resolve_file_values(val, dx_files_dir = 'input_files')

# ** _resolve_file_values

class FileResolver(object):
    """Abstract base class for recognizing and resolving values that point to files."""

    def __init__(self, input_files_dir):
        self.input_files_dir = input_files_dir
    
    def is_file_value(self, val): 
        """Test whether `val` represents a file that this FileResolver can handle"""
        raise NotImplementedError()

    def create_git_annex_symlink(self, file_val, git_file_path):
        """Given a value representing a file, create a git-annex link to the file under the relative path `git_file_path`.
        Returns the git-annex key representing the contents of `file_val`.
        """
        raise NotImplementedError()

class DxFileResolver(FileResolver):

    """FileResolver for DNAnexus links"""

    def is_file_value(self, val):
        """Test whether `val` represents a DNAnexus file"""
        return _is_mapping(val) and '$dnanexus_link' in val

    def _get_dx_id(self, file_val):
        """Given a `val` that represents a DNAnexus file (possibly indirectly), get the dx ID (file-xxxxxxx) for the file."""
        link = file_val['$dnanexus_link']
        if isinstance(link, collections.Mapping) and 'analysis' in link:
            # a link to the output of a particular stage of a particular DNAnexus analysis
            return self._get_dx_file_id(dxpy.DXAnalysis(link['analysis']).describe()['output'][link['stage']+'.'+link['field']])
        elif (_is_str(link) and link.startswith('file-')) or \
             (isinstance(link, collections.Mapping) and 'id' in link and \
              _is_str(link['id']) and link['id'].startswith('file-')):
            if _is_str(link) and link.startswith('file-'):
                return link
            else:
                return link['id']

    def create_git_annex_symlink(self, file_val, git_file_path):
        """Given a value representing a file, create a git-annex link to the file under the relative path `git_file_path`.
        Returns the git-annex key representing the contents of `file_val`.
        """
        dxid = self._get_dx_id(file_val)


        dx_descr = dxpy.describe(dxid)
        dx_uri = 'dx://' + dxid + '/' + descr['name']

        _run("git annex addurl --fast '{}' --file '{}' ".format(dx_uri, git_file_path))
        ga_key = _run_get_output("git annex lookupkey '{}'".format(git_file_path))

        def get_md5e_key(uri_key):
            uri_key_mdata = _run_get_json("git annex metadata --json --key='{}'".format(uri_key))
            if uri_key_mdata['fields']:
                md5e_keys = [k for k in uri_key_mdata['fields'].get('altKeys', ()) if k.startswith('MD5E-')]
                if md5e_keys:
                    assert len(md5e_keys) == 1
                    return md5e_keys[0]
            return None
        md5e_key = get_md5e_key(ga_key)
        if md5e_key:
            _run('git rm ' + git_file_path)
            _run("git annex fromkey '{}' '{}'".format(md5e_key, git_file_path))
            ga_key = md5e_key
        return ga_key
    # end: def create_git_annex_symlink(self, file_val, git_file_path):
# end: class DxFileResolver(FileResolver)

class GitAnnexTool(object):

    def get_metadata(self, key):
        """Return metadata associated with key, as a tuple of values.  If no metadata is associated
        with key, return empty tuple."""
        key_metadata = _run_get_json("git annex metadata --json --key='{}'".format(key))
        return tuple(key_metadata.get('fields', {}).get(key, ()))

class GitAnnexRepo(object):

    """A local git-annex repository, together with its configured cloud remotes.
    """

    def from_key(ga_key, repo_rel_path):
        """Set the file at `repo_rel_path` to point to the git-annex key `ga_key`."""
        pass

    def record_file_uri(file_uri):
        """Ensures that the given URI is recorded in this repo.

        Returns:
           A git-annex key for the contents of file_uri.  
        """
        pass

    pass

class GitAnnexRemote(object):

    def has_key(key):
        """Test whether the remote contains a given key"""
        raise NotImplementedError()

    def put_key(key):
        """Ensure the remote has content with given key."""
        raise NotImplementedError()

# ** _get_dx_val
def _get_dx_val(val, dx_files_dir):
    """Resolve a dx value: if it is a scalar, just return that;
    if it is a dx file, fetch the file, cache it locally, and return the path to the file."""

    #print('parsing val: ', val)
    util.file.mkdir_p(dx_files_dir)
    if isinstance(val, list):
        return [_get_dx_val(val=v, dx_files_dir=dx_files_dir) for v in val]
    if isinstance(val, collections.Mapping) and '$dnanexus_link' in val:
        link = val['$dnanexus_link']
        if isinstance(link, collections.Mapping) and 'analysis' in link:
            _log.debug('link is %s', link)
            return _get_dx_val(dxpy.DXAnalysis(link['analysis']).describe()['output'][link['stage']+'.'+link['field']],
                               dx_files_dir)
        elif (_is_str(link) and link.startswith('file-')) or (isinstance(link, collections.Mapping) and 'id' in link and _is_str(link['id']) and link['id'].startswith('file-')):
            if _is_str(link) and link.startswith('file-'):
                dxid = link
            else:
                dxid = link['id']
            descr = dxpy.describe(dxid)
            dx_file = os.path.join(dx_files_dir, dxid) + '-' + descr['name']
            file_size = int(descr['size'])

            # see if the file is cached in git-annex
            ga_mdata = _run_get_json('git annex metadata --json --key=WORM-s0-m0--dx-' + dxid)
            if ga_mdata['fields']:
                ga_key = ga_mdata['fields']['ga_key'][0]
                _run('git annex get --key', ga_key)
                _run('git annex fromkey', ga_key, dx_file)

            if not os.path.isfile(dx_file) or os.path.getsize(dx_file) != file_size:
                _log.debug('fetching %s to %s', dxid, dx_file)
                # TODO: check that there is enough free space (with some to spare)
                if os.path.isfile(dx_file):
                    os.unlink(dx_file)
                fs_info = os.statvfs(dx_files_dir)
                assert file_size < fs_info.f_bsize * fs_info.f_bavail
                dxpy.bindings.dxfile_functions.download_dxfile(dxid=dxid, filename=dx_file+'.fetching', show_progress=True)
                assert os.path.getsize(dx_file+'.fetching') == file_size
                os.rename(dx_file+'.fetching', dx_file)
                _log.debug('fetched %s to %s', dxid, dx_file)
                _log.debug('curdir is %s', os.getcwd())
                _run('git annex add', dx_file)
                ga_key = _run_get_output('git annex lookupkey ' + dx_file).strip()
                _run('git annex metadata -s', 'dxid+=' + dxid, dx_file)
                # record a mapping from the dxid to the git-annex key
                _run('git annex metadata', '--key=WORM-s0-m0-dx-'+dxid, '-s', 'ga_key='+ga_key)
                # register a URL that can be used to re-fetch this file from DNAnexus;
                # the URL is only valid if the 'run_dx_url_server' command is running.
                #_run('git annex registerurl ' + ga_key + ' ' + ' http://localhost:8080/dx/' + dxid)
        else:
            raise RuntimeError('Cannot parse dx link {}'.format(link))
        return dx_file
    else:
        return val
# end: def _get_dx_val(val, dx_files_dir = 'input_files')

def _construct_analysis_inputs(inputs):
    """"Given parsed specification of analysis inputs, construct for each analysis a dictionary of its specified inputs.
    Additionally, construct analysis labels which provide additional information about the analysis.
    """
    

# def _construct_run_inputs(workflow_name, inputs, analysis_dir):
#     """Construct the inputs for the given workflow.

#     Args:
#       workflow_name: name of the wdl workflow
#       inputs: ordered list of sources from which to take inputs.
#       analysis_dir: the analysis dir in which we are constructing the input spec
#       orig_cwd: the current working directory from which the top-level command was run

#     Returns:
#       a map from input name to a value.  values may include git_links to files under the analysis dir.
#     """
#     # workflow_defaults = { input_name: input_spec['default'] for input_name, input_spec in workflow_inputs_spec.items()
#     #                       if input_spec.get('default') is not None }

#     def _load_inputs(inputs_source):
#         # if inputs_source == '_workflow_defaults':
#         #     return workflow_defaults

#         key = None
#         if ':' in inputs_source:
#             inputs_source, key = inputs_source.split(':')
#         if not os.path.isabs(inputs_source):
#             inputs_source = os.path.join(orig_cwd, inputs_source)
#         if os.path.isdir(inputs_source) and os.path.isfile(os.path.join(inputs_source, 'metadata.json')):
#             inputs_source = os.path.join(inputs_source, 'metadata.json')
#         if key is None and os.path.basename(inputs_source) == 'metadata.json':
#             key = 'inputs'
#         inps = _json_loadf(inputs_source)
#         inps_orig = inps
#         if key is not None:
#             inps = inps[key]

#         # rename inputs, if needed
#         for k in set(inps):
#             for orig, repl in (input_name_subst or ()):
#                 if orig in k:
#                     _dict_rename_key(inps, k, k.replace(orig, repl))

#         # change the workflow name to ours, if needed
#         for k in set(inps):
#             if not k.startswith(workflow_name+'.'):
#                 assert k.startswith(inps_orig['workflowName']+'.')
#                 k_renamed=workflow_name+k[k.index('.'):]
#                 print('RENAMING INPUT KEY', k, k_renamed)
#                 _dict_rename_key(inps, k, k_renamed)

#         # change the paths in any git_links to be relative to the analysis dir
#         inputs_source_dir = os.path.dirname(inputs_source)
#         def reroute_git_link(val):
#             if not _maps(val, '$git_link'):
#                 return val
#             return {'$git_link': os.path.relpath(os.path.join(inputs_source_dir, val['$git_link']),
#                                                  analysis_dir)}
#         return _transform_json_data(inps, reroute_git_link)

#     inp_srcs = list(map(_load_inputs, inputs_sources))
#     run_inputs = _ord_dict()
#     for inp_src in inp_srcs:
#         run_puts.update(inp_src)

#     print('CONSTRUCTED', run_inputs)
#     return run_inputs
# end: def _construct_run_inputs(workflow_name, inputs_sources, input_name_subst, analysis_dir)

# ** deleted code

# def _resolve_links_in_dx_analysis(dx_analysis_id, analysis_dir):
#     analysis_descr = _dx_describe(dx_analysis_id)
#     del analysis_descr['originalInput']
#     methods = [functools.partial(_resolve_link_dx, dx_analysis_id=dx_analysis_id)]
#     resolved = _resolve_links_in_json_data(val=analysis_descr, rel_to_dir=analysis_dir, methods=methods)
#     _write_json(os.path.join(analysis_dir, 'dx_resolved.json'), **resolved)

# def _resolve_link_dx_old(val, dx_analysis_id):
#     """Resolve a dx value: if it is a scalar, just return that;
#     if it is a dx link, resolve the link to a file-xxxx identifier."""

#     print('parsing val: ', val)
#     if not (_is_mapping(val) and '$dnanexus_link' in val):
#         return None

#     recurs = functools.partial(_resolve_dx_ids_in_val, analysis_descr=analysis_descr)
#     if isinstance(val, list):
#         return [recurs(val=v, git_path = os.path.join(git_path, str(i))) for i, v in enumerate(val)]
#     if isinstance(val, collections.Mapping):
#         if  '$dnanexus_link' not in val:
#             return collections.OrderedDict([(k, recurs(val=v, git_path=os.path.join(git_path, util.file.string_to_file_name(k))))
#                                             for k, v in val.items()])
#         link = val['$dnanexus_link']
#         if isinstance(link, collections.Mapping) and 'stage' in link and ('field' in link or 'outputField' in link):
#             print('link is ', link)
#             linked_analysis_descr = dxpy.DXAnalysis(link['analysis']).describe() if 'analysis' in link else analysis_descr
#             linked_field = link['field'] if 'field' in link else link['outputField']
#             return recurs(val=linked_analysis_descr['output'][link['stage']+'.'+linked_field],
#                           analysis_descr=analysis_descr, git_path=git_path)
#         elif (_is_str(link) and link.startswith('file-')) or (isinstance(link, collections.Mapping) and 'id' in link and _is_str(link['id']) and link['id'].startswith('file-')):
#             if _is_str(link) and link.startswith('file-'):
#                 dxid = link
#             else:
#                 dxid = link['id']




#     def _resolve_dx_ids_in_val(val, analysis_descr, git_path):
#         """Resolve a dx value: if it is a scalar, just return that;
#         if it is a dx link, resolve the link to a file-xxxx identifier."""

#         print('parsing val: ', val)

#         recurs = functools.partial(_resolve_dx_ids_in_val, analysis_descr=analysis_descr)
#         if isinstance(val, list):
#             return [recurs(val=v, git_path = os.path.join(git_path, str(i))) for i, v in enumerate(val)]
#         if isinstance(val, collections.Mapping):
#             if  '$dnanexus_link' not in val:
#                 return collections.OrderedDict([(k, recurs(val=v, git_path=os.path.join(git_path, util.file.string_to_file_name(k))))
#                                                 for k, v in val.items()])
#             link = val['$dnanexus_link']
#             if isinstance(link, collections.Mapping) and 'stage' in link and ('field' in link or 'outputField' in link):
#                 print('link is ', link)
#                 linked_analysis_descr = dxpy.DXAnalysis(link['analysis']).describe() if 'analysis' in link else analysis_descr
#                 linked_field = link['field'] if 'field' in link else link['outputField']
#                 return recurs(val=linked_analysis_descr['output'][link['stage']+'.'+linked_field],
#                               analysis_descr=analysis_descr, git_path=git_path)
#             elif (_is_str(link) and link.startswith('file-')) or (isinstance(link, collections.Mapping) and 'id' in link and _is_str(link['id']) and link['id'].startswith('file-')):
#                 if _is_str(link) and link.startswith('file-'):
#                     dxid = link
#                 else:
#                     dxid = link['id']
#                 util.file.mkdir_p(git_path)
#                 git_file_path = import_from_url(url='dx://' + dxid, git_file_path = git_path)
#                 return _ord_dict(('$git_relpath', os.path.relpath(git_file_path, )) # return { $git_link: }
#                 descr = dxpy.describe(dxid)
#                 dx_file = os.path.join(dx_files_dir, dxid) + '-' + descr['name']
#                 file_size = int(descr['size'])

#                 # see if the file is cached in git-annex
#                 ga_mdata = _run_get_json('git annex metadata --json --key=WORM-s0-m0--dx-' + dxid)
#                 if ga_mdata['fields']:
#                     ga_key = ga_mdata['fields']['ga_key'][0]
#                     _run('git annex get --key', ga_key)
#                     _run('git annex fromkey', ga_key, dx_file)

#                 if not os.path.isfile(dx_file) or os.path.getsize(dx_file) != file_size:
#                     print('fetching', dxid, 'to', dx_file)
#                     # TODO: check that there is enough free space (with some to spare)
#                     if os.path.isfile(dx_file):
#                         os.unlink(dx_file)
#                     fs_info = os.statvfs(dx_files_dir)
#                     assert file_size < fs_info.f_bsize * fs_info.f_bavail
#                     dxpy.bindings.dxfile_functions.download_dxfile(dxid=dxid, filename=dx_file+'.fetching', show_progress=True)
#                     assert os.path.getsize(dx_file+'.fetching') == file_size
#                     os.rename(dx_file+'.fetching', dx_file)
#                     print('fetched', dxid, 'to', dx_file)
#                     print('curdir is ', os.getcwd())
#                     _run('git annex add', dx_file)
#                     ga_key = _run_get_output('git annex lookupkey ' + dx_file).strip()
#                     _run('git annex metadata -s', 'dxid+=' + dxid, dx_file)
#                     # record a mapping from the dxid to the git-annex key
#                     _run('git annex metadata', '--key=WORM-s0-m0-dx-'+dxid, '-s', 'ga_key='+ga_key)
#                     # register a URL that can be used to re-fetch this file from DNAnexus;
#                     # the URL is only valid if the 'run_dx_url_server' command is running.
#                     #_run('git annex registerurl ' + ga_key + ' ' + ' http://localhost:8080/dx/' + dxid)
#             else:
#                 raise RuntimeError('Cannot parse dx link {}'.format(link))
#             return dx_file
#         else:
#             return val
# # end: def _resolve_dx_val(val, analysis_descr, dx_files_dir = 'input_files')




# def _record_dx_metadata(val, analysis_dir, root_dir):
#     """If 'val' is a filename, return a dict representing the file and some metadata about it;
#     otherwise, return val as-is."""
#     if isinstance(val, list): return [_record_file_metadata(v, analysis_dir, root_dir) for v in val]
#     if _is_mapping(val): return collections.OrderedDict([(k, _record_file_metadata(v, analysis_dir, root_dir))
#                                                           for k, v in val.items()])
#     if not (_is_str(val) and (os.path.isfile(val) or os.path.isdir(val))): return val
#     file_info = collections.OrderedDict([('_is_file' if os.path.isfile(val) else '_is_dir', True)])
#     assert val.startswith(analysis_dir) or val.startswith(root_dir)
#     if val.startswith(analysis_dir):
#         relpath = os.path.relpath(val, analysis_dir)
#         abspath = os.path.join(analysis_dir, relpath)
#     else:
#         cromwell_executions_dir = os.path.dirname(os.path.dirname(root_dir))
#         relpath = os.path.relpath(val, cromwell_executions_dir)
#         abspath = os.path.join(analysis_dir, 'output',
#                                'call_logs' if os.path.basename(val) in ('stdout', 'stderr') else 'outputs', relpath)
#         if os.path.isfile(val) and not os.path.isfile(abspath):
#             print('LINKING {} to {}'.format(val, abspath))
#             util.file.mkdir_p(os.path.dirname(abspath))
#             shutil.copy(val, abspath)
#         if os.path.isdir(val):
#             util.file.mkdir_p(abspath)

#     assert os.path.isabs(abspath) and abspath.startswith(analysis_dir), \
#         'bad abspath: {} analysis_dir: {}'.format(abspath, analysis_dir)
#     relpath = os.path.relpath(abspath, analysis_dir)
#     assert not os.path.isabs(relpath), 'should be relative: {}'.format(relpath)
#     assert os.path.isfile(abspath) or os.path.isdir(abspath), 'not file or dir: {}'.format(abspath)
#     assert os.path.isdir(abspath) or os.path.getsize(abspath) == os.path.getsize(val)
#     file_info['relpath'] = relpath
#     if os.path.isfile(abspath):
#         file_info['size'] = os.path.getsize(abspath)
#         file_info['md5'] = _run_get_output('md5sum ' + abspath).strip().split()[0]
#     return file_info

def _resolve_link_gs(val, git_file_dir, fast=True):
    if not (_is_str(val) and val.startswith('gs://') and _gs_stat(val)): return val
    util.file.mkdir_p(git_file_dir)
    val_resolved = {'$git_link': import_from_url(url=val, git_file_path=git_file_dir, fast=fast)}
    return val_resolved

def _resolve_link_local_path(val, git_file_dir):
    if not (_is_str(val) and os.path.lexists(val)):
        return val
    util.file.mkdir_p(git_file_dir)
    fname = os.path.join(git_file_dir, os.path.basename(val))
    if os.path.lexists(fname):
        os.unlink(fname)
        
    link_into_annex, link_target_in_annex = tools.git_annex.GitAnnexTool()._get_link_into_annex(val)
    if link_target_in_annex:
        os.symlink(link_into_annex, fname)
    elif os.path.isfile(val):
        shutil.copyfile(val, fname)
        assert os.path.isfile(fname)
        _git_annex_add(fname)
    else:
        return val
    _log.debug('GIT LINK TO %s', fname)
    return {'$git_link': fname}
    
    # for now can just copy (or hardlink? careful) it here and git-annex-add
    # later can see if can use known md5 from cromwell's caching calculations

    # what if it's a dir?  
def _resolve_link_dx2(val, git_file_dir, dx_analysis_id, _cache=None, git_annex_tool=None):
    val_resolved = _resolve_dx_link_to_dx_file_id_or_value(val=val, dx_analysis_id=dx_analysis_id)
    if _maps(val_resolved, '$dnanexus_link'):
        dx_file_id = val_resolved['$dnanexus_link']
        assert dx_file_id.startswith('file-')
        if _maps(_cache, dx_file_id):
            val_resolved = copy.deepcopy(_cache[dx_file_id])
            #_log.debug('RESOLVED %s TO %s FROM CACHE', val, val_resolved)
        else:
            util.file.mkdir_p(git_file_dir)
            # git_annex_tool.  *** CHANGE HERE
            val_resolved = {'$git_link': import_from_url(url='dx://' + dx_file_id, git_file_path=git_file_dir, fast=False)}
            #_log.debug('RESOLVED %s TO %s', val, val_resolved)
            if _cache is not None:
                _cache[dx_file_id] = val_resolved
    return val_resolved

def _resolve_link_dx(val, git_file_dir, dx_analysis_id, _cache=None):
    val_resolved = _resolve_dx_link_to_dx_file_id_or_value(val=val, dx_analysis_id=dx_analysis_id)
    if _maps(val_resolved, '$dnanexus_link'):
        dx_file_id = val_resolved['$dnanexus_link']
        assert dx_file_id.startswith('file-')
        if _maps(_cache, dx_file_id):
            val_resolved = copy.deepcopy(_cache[dx_file_id])
            #_log.debug('RESOLVED %s TO %s FROM CACHE', val, val_resolved)
        else:
            util.file.mkdir_p(git_file_dir)
            val_resolved = {'$git_link': import_from_url(url='dx://' + dx_file_id, git_file_path=git_file_dir, fast=False)}
            #_log.debug('RESOLVED %s TO %s', val, val_resolved)
            if _cache is not None:
                _cache[dx_file_id] = val_resolved
    return val_resolved

# * RedirectToCloudObject

class RedirectToCloudObject(SimpleHTTPServer.SimpleHTTPRequestHandler):

   def _get_dx(self):
       try:
           file_info = _run_get_json('dx', 'describe', '--json', self._get_dx_id())
       except subprocess.CalledProcessError as e:
           _log.info('CalledProcessError: %s', e)
           raise
       _log.debug('got file_info %s', file_info)
       assert file_info['id'] == self._get_dx_id()

       self.send_response(307)
       if False and 'media' in file_info:
           self.send_header('Content-Type', file_info['media'])
       if False and 'size' in file_info:
           self.send_header('Content-Length', str(file_info['size']))
       #self.send_header('Content-Disposition', 'attachment; filename="{}"'.format(file_info['name']))

       new_path = _run_get_output('dx make_download_url --duration 2h ' + self._get_dx_id())
       _log.debug('new_path=%s', new_path)
       self.send_header('Location', new_path)
       self.end_headers()


   def _get_gs(self):
       _log.info('IN_gs_gs')
       try:
           _log.info('gs_url is %s', self._get_gs_uri())
           signed_url = _run_get_output('gsutil signurl -d 10m ' + self.server.gs_key + ' ' + self._get_gs_uri()).strip().split()[-1]
           _log.info('signed url is %s', signed_url)
       except subprocess.CalledProcessError as e:
           _log.error('CalledProcessError: %s', e)
           raise

       self.send_response(307)
       if False and 'media' in file_info:
           self.send_header('Content-Type', file_info['media'])
       if False and 'size' in file_info:
           self.send_header('Content-Length', str(file_info['size']))
       #self.send_header('Content-Disposition', 'attachment; filename="{}"'.format(file_info['name']))

       self.send_header('Location', signed_url)
       self.end_headers()

   def do_GET(self):
       try:
           _log.info('path is %s', self.path)
           _log.info('is gs path? %s', self.path.startswith('/gs/'))
           if not any(self.path.startswith('/' + p) for p in self.server.valid_paths):
               raise IOError('not a valid path')
           if self.path.startswith('/dx/'):
               self._get_dx()
           elif self.path.startswith('/gs/') and self.server.gs_key:
               _log.info('CALLING get_gs')
               self._get_gs()
               _log.info('RETURNED FROM get_gs')
           else:
               _log.warning('UNKNOWN PATH!! %s', self.path)
               raise IOError('Unknown URL')
       except (IOError, subprocess.CalledProcessError):
           self.send_error(404, 'file not found')

   def do_HEAD(self):
       self.do_GET()

   def _get_dx_id(self):
       return self.path[len('/dx/'):]

   def _get_gs_uri(self):
       return 'gs://' + self.path[len('/gs/'):]

def run_cloud_object_url_server(port, gs_key, valid_paths):
    """Start a webserver that will redirect URL requests of the form http://localhost/dx/file-xxxxxx to DNAnexus,
    requests of the form http://localhost/gs/gsbucket/gsobjectpath to Google Cloud Storage, and requests
    of the form http://localhost/s3/s3bucket/gsobjectpath to AWS S3.
    This gives each cloud object a stable http URL, permitting such files to be added to git-annex with the addurl command.
    """

    assert valid_paths, 'Some valid paths must be specified'
    server = None
    try:           
        SocketServer.TCPServer.allow_reuse_address = True
        server = SocketServer.TCPServer(("", port), RedirectToCloudObject)
        server.gs_key = gs_key
        server.valid_paths = valid_paths
        server.serve_forever()
    except Exception as e:
        if server is None:
            _log.error('No server!')
        else:
            _log.info('calling shutdown...')
            server.shutdown()
            _log.info('calling server_close...')
            server.server_close()
            _log.error('re-raising exception %s', e)
            raise

def parser_run_cloud_object_url_server(parser=argparse.ArgumentParser()):
    parser.add_argument('--port', default=8080, help='Port on which to run the webserver')
    parser.add_argument('--gsKey', dest='gs_key', help='Key for signing gs urls')
    parser.add_argument('--validPaths', dest='valid_paths', required=True, help='paths that may be accessed',
                        nargs='+')
    util.cmd.attach_main(parser, run_cloud_object_url_server, split_args=True)
    return parser

__commands__.append(('run_cloud_object_url_server', parser_run_cloud_object_url_server))

########################################################################################################################
