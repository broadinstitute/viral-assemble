"""Code for comparing assemblies"""

import os
import asmutil

class AsmCmp(object):
    """Code for comparing assemblies"""

    def __init__( self, cfg ):
        self.cfg = cfg

    def put_stat( self, readset, assembler, metric, value ):
        fname = self.__get_stat_fname( readset, assembler, metric )
        asmutil.ensure_dir_exists( os.path.dirname( fname ) )
        asmutil.dump_file( fname, value )

    def get_stat( self, readset, assembler, metric ):
        FN = self.__get_stat_fname( readset, assembler, metric )
        return float( asmutil.slurp_file( FN ) ) if  os.path.exists( FN ) else None

    def __get_stat_fname( self, readset, assembler, metric ):
        return os.path.join( self.cfg['statdir'], readset, assembler, metric + '.dat' )

# end: class AsmCmp

def test_AsmCmp():
    statdir = '/broad/hptmp/ilya/asmstats'
    asmutil.ensure_dir_exists( statdir )
    cfg = dict( statdir = statdir )
    asmcmp = AsmCmp( cfg = cfg )
    asmcmp.put_stat( 'READS1', 'NOVALIGN_1', 'NA50', 1000 )
    print( asmcmp.get_stat( 'READS1', 'NOVALIGN_1', 'NA50' ) )

if __name__ == '__main__':
    test_AsmCmp()


        
