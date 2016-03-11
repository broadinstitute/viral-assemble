import logging
import tools
import util.file
import os
import os.path
import subprocess
from collections import OrderedDict
#import pysam

TOOL_NAME = 'sambamba'
TOOL_VERSION = '0.5.9'
TOOL_URL = 'http://sourceforge.net/projects/samtools/files/samtools/{ver}/samtools-{ver}.tar.bz2'.format(
    ver=TOOL_VERSION)

TOOL_URL = 'https://github.com/lomereiter/sambamba/releases/download/v{ver}/sambamba_v{ver}_{ostype}.tar.bz2'

log = logging.getLogger(__name__)


class SambambaTool(tools.Tool):

    def __init__(self, install_methods=None):
        if install_methods is None:
            install_methods = []
            install_methods.append(tools.CondaPackage(TOOL_NAME, version=TOOL_VERSION))
            install_methods.append(
                tools.DownloadPackage(
                    TOOL_URL.format(ver=TOOL_VERSION,
                                    ostype=get_sambamba_os()),
                    'sambamba_v{}'.format(TOOL_VERSION),
                    verifycmd='{}/sambamba_v{} view &> /dev/null'.format(util.file.get_build_path(), TOOL_VERSION)
                )
            )
        tools.Tool.__init__(self, install_methods=install_methods)

    def version(self):
        return TOOL_VERSION

    def execute(self, command, args, stdin=None, stdout=None, stderr=None):    # pylint: disable=W0221
        tool_cmd = [self.install_and_get_path(), command] + args
        log.debug(' '.join(tool_cmd))
        if stdin:
            stdin = open(stdin, 'r')
        if stdout:
            stdout = open(stdout, 'w')
        if stderr:
            stderr = open(stderr, 'w')
        subprocess.check_call(tool_cmd, stdin=stdin, stdout=stdout, stderr=stderr)
        if stdin:
            stdin.close()
        if stdout:
            stdout.close()
        if stderr:
            stderr.close()

    def view(self, args, in_file, out_file, regions=None, threads=None):
        regions = regions or []
        opts = []

        if threads:
            opts = opts + ['-t', threads]

        self.execute('view', args + opts + ['-o', out_file, in_file] + regions)

    def count(self, in_bam, opts=None, regions=None, threads=None):
        opts = opts or []
        regions = regions or []

        if threads:
            opts = opts + ['-t', threads]

        cmd = [self.install_and_get_path(), 'view', '-c'] + opts + [in_bam] + regions
        # return int(pysam.view(*cmd)[0].strip())
        return int(subprocess.check_output(cmd).strip())

    def index(self, in_bam, threads=None):
        opts = []

        if threads:
            opts = opts + ['-t', threads]

        self.execute('index', opts + [in_bam])

    def merge(self, in_files, out_file, options=None):
        "Merge a list of in_files to create out_file."
        options = options or []

        self.execute('merge', options + [out_file] + in_files)

    def sort(self, in_file, out_file, args=None):
        args = args or []

        self.execute('sort', args + ['-o', out_file, in_file])

    def flagstat(self, in_bam):
        self.execute('flagstat', [in_bam])

    def slice(self, in_file, out_file, region=None, threads=None):
        opts = []
        if threads:
            opts = opts + ['-t', threads]

        self.execute('slice', opts + ['-o', out_file, in_file] + region)

    def markdup(self, in_bam, out_bam):
        opts = ['-r']
        self.execute('sort', opts + [in_bam, out_bam])

    def depth(self, in_file):
        raise NotImplementedError("sambamba depth not yet implemented")

    def mpileup(self, in_bam, outPileup, opts=None, threads=None):
        # requires samtools and bcftools
        raise NotImplementedError("sambamba mpileup not yet implemented")

        opts = opts or []

        if threads:
            opts = opts + ['-t', threads]

        self.execute('mpileup', opts + [in_bam], stdout=outPileup, stderr='/dev/null')    # Suppress info messages

    def dumpHeader(self, in_bam, outHeader):
        if in_bam.endswith('.bam'):
            opts = ['-H']
        elif in_bam.endswith('.sam'):
            opts = ['-H', '-S']
        #header = pysam.view(*opts)
        self.view(opts, in_bam, outHeader)

    def getHeader(self, in_bam):
        ''' fetch BAM header as a list of tuples (already split on tabs) '''
        tmpf = util.file.mkstempfname('.txt')
        self.dumpHeader(in_bam, tmpf)
        with open(tmpf, 'rb') as inf:
            header = list(line.decode("latin-1").rstrip('\n').split('\t') for line in inf)
        os.unlink(tmpf)
        return header

    def getReadGroups(self, in_bam):
        ''' fetch all read groups from the BAM header as an OrderedDict of
            RG ID -> RG dict.  The RG dict is a mapping of read group keyword
            (like ID, DT, PU, LB, SM, CN, PL, etc) to value.  ID is included
            and not stripped out. ID is required for all read groups.
            Resulting keys are in same order as @RG lines in bam file.
        '''
        rgs = [
            dict(x.split(':', 1) for x in row[1:]) for row in self.getHeader(in_bam) if len(row) > 0 and row[0] == '@RG'
        ]
        return OrderedDict((rg['ID'], rg) for rg in rgs)

    # left over from samtools...

    def faidx(self, inFasta, overwrite=False):
        raise NotImplementedError("not yet implemented")

    def reheader(self, in_bam, headerFile, out_bam):
        raise NotImplementedError("not yet implemented")

    def removeDoublyMappedReads(self, in_bam, out_bam):
        raise NotImplementedError("not yet implemented")


def get_sambamba_os():
    uname = os.uname()
    if uname[0] == "Darwin":
        return "osx"
    if uname[0] == "Linux":
        return "linux"
