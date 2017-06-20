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
        install_methods = [tools.CondaPackage(
            TOOL_NAME,
            version=TOOL_VERSION,
            executable='gatk-launch')]
        tools.Tool.__init__(self, install_methods=install_methods)

    def execute(self, command, gatkOptions=None, JVMmemory=None):    # pylint: disable=W0221
        gatkOptions = gatkOptions or []

        if not JVMmemory:
            JVMmemory = self.jvmMemDefault

        # the conda version wraps the jar file with a shell script
        tool_cmd = [
            self.install_and_get_path(),
            '--sparkRunner', 'LOCAL',
            '--javaOptions', '-Xmx%s -Djava.io.tmpdir=%s' % (JVMmemory, tempfile.gettempdir()),
            command
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
         self.tool_version = TOOL_VERSION

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
