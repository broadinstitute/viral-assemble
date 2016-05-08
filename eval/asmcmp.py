"""Code for comparing assemblies"""

import os
import asmutil

class AsmCmp(object):
    """Code for comparing assemblies"""

    def __init__( self, cfg ):
        self.cfg = cfg

    def put_stat( self, sample, asm, stat, value ):
        fname = self.__get_stat_fname( sample, asm, stat )
        asmutil.ensure_dir_exists( os.path.dirname( fname ) )
        asmutil.dump_file( fname, value )

    def get_stat( self, sample, asm, stat ):
        FN = self.__get_stat_fname( sample, asm, stat )
        return float( asmutil.slurp_file( FN ) ) if  os.path.exists( FN ) else None

    def __get_stat_fname( self, sample, asm, stat ):
        """Return the name of the file where the given statistic is stored"""
        return os.path.join( self.cfg['statdir'], sample, asm, stat + '.dat' )

    def find_largest_changes( self, asm1, asm2, samples, stats, out_dir ):
        """For a pair of assemblers, find the largest differences between their performance."""
        
        for stat in stats:
            deltas = []
            for sample in samples:
                stat_asm1 = self.get_stat( sample, asm1, stat )
                stat_asm2 = self.get_stat( sample, asm2, stat )
                if stat_asm1 is not None and stat_asm2 is not None:
                    deltas.append( stat_asm1 - stat_asm2 )

            import matplotlib as mp
            mp.use('agg')
            import matplotlib.pyplot as pp

            pp.figure()
            pp.xlabel( 'delta ' + stat )
            pp.title( asm1 + ' vs. ' + asm2 )
            pp.hist(deltas)
            asmutil.ensure_dir_exists( out_dir )
            pp.savefig( os.path.join( out_dir, stat + '.png' ) )

# end: class AsmCmp

def test_AsmCmp():
    statdir = '/broad/hptmp/ilya/asmstats'
    asmutil.ensure_dir_exists( statdir )
    cfg = dict( statdir = statdir )
    asmcmp = AsmCmp( cfg = cfg )
    asmcmp.put_stat( 'READS1', 'NOVALIGN_1', 'NA50', 1000 )
    asmcmp.put_stat( 'READS2', 'NOVALIGN_1', 'NA50', 2000 )
    asmcmp.put_stat( 'READS1', 'NOVALIGN_2', 'NA50', 500 )
    asmcmp.put_stat( 'READS2', 'NOVALIGN_2', 'NA50', 700 )
    #print( asmcmp.get_stat( 'READS1', 'NOVALIGN_1', 'NA50' ) )
    asmcmp.find_largest_changes( asm1 = 'NOVALIGN_1', asm2 = 'NOVALIGN_2', 
                                 samples = ( 'READS1', 'READS2' ),
                                 stats = [ 'NA50' ], out_dir = os.path.join( statdir, 'plots' ) )

if __name__ == '__main__':
    test_AsmCmp()


        
