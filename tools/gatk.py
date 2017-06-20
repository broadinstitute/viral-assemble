'''
    GATK4 genotyping toolkit from the Broad Institute
'''

import tools
import tools.picard
import tools.samtools
import util.file
import util.misc

import logging
import os
import os.path
import subprocess
import tempfile

_log = logging.getLogger(__name__)


TOOL_NAME = 'gatk4'
TOOL_VERSION = '4.0a1.2.7.2'

class GATKTool(tools.Tool):
    jvmMemDefault = '2g'

    def __init__(self):
        self.tool_version = None
        install_methods = [tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION)]
        tools.Tool.__init__(self, install_methods=install_methods)

    def execute(self, command, gatkOptions=None, JVMmemory=None):    # pylint: disable=W0221
        gatkOptions = gatkOptions or []

        if not JVMmemory:
            JVMmemory = self.jvmMemDefault

        # the conda version wraps the jar file with a shell script
        if self.install_and_get_path().endswith(".jar"):
            tool_cmd = [
                'java', '-Xmx' + JVMmemory, '-Djava.io.tmpdir=' + tempfile.gettempdir(), '-jar', self.install_and_get_path(),
                '-T', command
            ] + list(map(str, gatkOptions))
        else:
            tool_cmd = [
                self.install_and_get_path(), '-Xmx' + JVMmemory, '-Djava.io.tmpdir=' + tempfile.gettempdir(), '-T', command
            ] + list(map(str, gatkOptions))

        _log.debug(' '.join(tool_cmd))
        subprocess.check_call(tool_cmd)

    @staticmethod
    def dict_to_gatk_opts(options):
        return ["%s=%s" % (k, v) for k, v in options.items()]

    def version(self):
        if self.tool_version is None:
            self._get_tool_version()
        return self.tool_version

    def _get_tool_version(self):
        if self.install_and_get_path().endswith(".jar"):
            cmd = [
                'java', '-Djava.io.tmpdir=' + tempfile.gettempdir(), '-jar', self.install_and_get_path(),
                '--version'
            ]
        else:
            cmd = [
                self.install_and_get_path(), '--version'
            ]

        self.tool_version = util.misc.run_and_print(cmd, buffered=False, silent=True).stdout.decode("utf-8").strip()

    def ug(self, inBam, refFasta, outVcf, options=None, JVMmemory=None, threads=1):
        options = options or ["--min_base_quality_score", 15, "-ploidy", 4]

        if int(threads) < 1:
            threads = 1
        opts = [
            '-I',
            inBam,
            '-R',
            refFasta,
            '-o',
            outVcf,
            '-glm',
            'BOTH',
            '--baq',
            'OFF',
            '--useOriginalQualities',
            '-out_mode',
            'EMIT_ALL_SITES',
            '-dt',
            'NONE',
            '--num_threads',
            threads,
            '-stand_call_conf',
            0,
            '-stand_emit_conf',
            0,
            '-A',
            'AlleleBalance',
        ]
        self.execute('UnifiedGenotyper', opts + options, JVMmemory=JVMmemory)

    def local_realign(self, inBam, refFasta, outBam, JVMmemory=None, threads=1):
        intervals = util.file.mkstempfname('.intervals')
        opts = ['-I', inBam, '-R', refFasta, '-o', intervals]
        _log.debug("Running local realign with %s threads", threads)
        self.execute('RealignerTargetCreator', opts, JVMmemory=JVMmemory)
        opts = ['-I',
                inBam,
                '-R',
                refFasta,
                '-targetIntervals',
                intervals,
                '-o',
                outBam,    #'--num_threads', threads,
               ]
        self.execute('IndelRealigner', opts, JVMmemory=JVMmemory)
        os.unlink(intervals)
