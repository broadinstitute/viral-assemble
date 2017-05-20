#!/usr/bin/env python
''' Commands for managing viral-ngs projects.  

See also https://github.com/broadinstitute/viral-ngs-deploy/tree/master/easy-deploy-script .
'''

__author__ = "ilya@broadinstitute.org"
__commands__ = []

import argparse, logging, os, yaml
import util.cmd, util.file, util.misc

log = logging.getLogger(__name__)

# =======================

def add_project_method(method_id, method_bin_dir, parent_method="default", 
                       inherit_results='data/source:.bam,data/depletion:.cleaned.bam,data/depletion:.taxfilt.bam,data/per_sample:cleaned.bam,data/per_sample:.taxfilt.bam'):
    """Set up directories for analyzing the current project by a different method.
    Results for method X are stored under methods/X/... .  Each method is defined by a config file.

    Args:
       method_id: an ID for this analysis method
       method_bin_dir: path to the viral-ngs code for this method
       parent_method: method from which we inherit the config file, and possibly some samples
       inherit_results: if given, hardlink the output for `inherit_samples` for these result steps.
    """

    assert method_id!='default'

    with open("config.yaml") as cfg_f:
        cfg = yaml.safe_load(f)

    methods_dir=cfg.get("methods_dir","methods")
    util.file.mkdir_p(os.path.join(methods_dir, method_id))
    if not os.path.exists(os.path.join(methods_dir, 'default')):
        os.symlink(os.path.realpath('.'), os.path.join(methods_dir, 'default'))

    bin_new=os.path.join(methods_dir, method_id, "bin")
    if os.path.exists(bin_new): 
        log.info('removing previous bin link {}'.format(bin_new))
        os.unlink(bin_new)
    os.symlink(method_bin_dir, bin_new)
    log.info('linked bin dir {} to {}'.format(method_bin_dir, os.path.join(methods_dir, "bin")))

    for subdir in cfg["subdirs"].values():
        util.file.mkdir_p(os.path.join(methods_dir, method_id, cfg["data_dir"], subdir))
        util.file.mkdir_p(os.path.join(methods_dir, method_id, cfg["tmp_dir"], subdir))
    util.file.mkdir_p(os.path.join(methods_dir, method_id, cfg["log_dir"]))
    util.file.mkdir_p(os.path.join(methods_dir, method_id, cfg["reports_dir"]))

    inherit_samples = util.file.slurp_file(cfg['samples_assembly']).strip().split()

    for sample in inherit_samples:
        for cfg_path_ext in inherit_results.split(','):
            cfg_path, ext = cfg_path_ext.split(':')
            cfg_path_dirs = [ cfg[d] for d in cfg_path.split('/') ]

            fn_old = os.path.join(*([methods_dir, parent_method] + cfg_path_dirs + [sample+ext]))
            fn_new = os.path.join(*([methods_dir, method_id ] + cfg_path_dirs + [sample+ext]))

            if os.path.isfile(fn_old) and not os.path.isfile(fn_new):
                log.info('linking {} to {}'.format(fn_old, fn_new))
                util.file.mkdir_p(os.path.dirname(fn_new))
                os.link(fn_old, fn_new)

    cfg_old = os.path.join(methods_dir, parent_method, 'config.yaml')
    cfg_new = os.path.join(methods_dir, method_id, 'config.yaml')
    log.info('cfg_old={} cfg_old_exists={}'.format(cfg_old, os.path.isfile(cfg_old)))
    log.info('cfg_new={} cfg_new_exists={}'.format(cfg_new, os.path.isfile(cfg_new)))
    if os.path.isfile(cfg_old) and not os.path.isfile(cfg_new):
        with open(cfg_old) as f_old, open(cfg_new, 'wt') as f_new:
            for line in f_old:
                if line.startswith('method_id:'):
                    line = 'method_id: "{}"'.format(method_id)
                for d in 'data tmp log reports bin'.split():
                    if line.startswith(d+'_dir'):
                        line = d+'_dir: "{}"'.format(os.path.join(methods_dir, method_id, cfg[d+'_dir']))
                f_new.write(line)
        log.info('rewrote config {} to {}'.format(cfg_old, cfg_new))
    else:
        log.info('not rewriting config')
        

def parser_add_project_method(parser=argparse.ArgumentParser()):
    parser.add_argument('method_id', help='Assembly method ID.')
    parser.add_argument('method_bin_dir', help='Location of viral-ngs code implementing the method.')
    parser.add_argument('--parent_method', help='Parent method from which to inherit config and partial results.', 
                        default='default')
    parser.add_argument('--inherit_results', default='data/source:.bam,data/depletion:.cleaned.bam,data/depletion:.taxfilt.bam,data/per_sample:cleaned.bam,data/per_sample:.taxfilt.bam', help='Results to inherit from parent')

    util.cmd.common_args(parser, (('loglevel', None), ('version', None)))
    util.cmd.attach_main(parser, add_project_method, split_args=True)
    return parser

__commands__.append(('add_project_method', parser_add_project_method))

# =======================


def full_parser():
    return util.cmd.make_parser(__commands__, __doc__)


if __name__ == '__main__':
    util.cmd.main_argparse(__commands__, __doc__)
