"""Code for comparing assemblies"""

from __future__ import division
import os, pandas, logging
import asmutil

class AsmCmp(object):
    """Code for comparing assemblies"""

    def __init__( self, cfg ):
        self.cfg = cfg

    def put_stat( self, readset, asm, stat, value ):
        fname = self.__get_stat_fname( readset, asm, stat )
        asmutil.ensure_dir_exists( os.path.dirname( fname ) )
        asmutil.dump_file( fname, value )

    def get_stat( self, readset, asm, stat ):
        FN = self.__get_stat_fname( readset, asm, stat )
        return float( asmutil.slurp_file( FN ) ) if  os.path.exists( FN ) else None

    def __get_stat_fname( self, readset, asm, stat ):
        """Return the name of the file where the given statistic is stored"""
        return os.path.join( self.cfg['statdir'], asmutil.make_alpha_num( readset ), 
                             asmutil.make_alpha_num( asm ), asmutil.make_alpha_num( stat ) + '.dat' )

    def find_largest_changes( self, asm1, asm2, readsets, stats, out_dir ):
        """For a pair of assemblers, find the largest differences between their performance."""
        
        for stat in stats:
            logging.info('find_largest_changes: stat=%s', stat)
            deltas = []
            for readset in readsets:
                stat_asm1 = self.get_stat( readset, asm1, stat )
                stat_asm2 = self.get_stat( readset, asm2, stat )
                if stat_asm1 is not None and stat_asm2 is not None:
                    diff = stat_asm1 - stat_asm2
                    if diff != 0:
                        logging.info( 'readset=%s stat=%s diff=%f', readset, stat, diff )
                    deltas.append( diff  )
            if not deltas: 
                logging.warn('no stats for ' + stat)
                continue
            import matplotlib as mp
            mp.use('agg')
            import matplotlib.pyplot as pp

            fig=pp.figure()
            pp.xlabel( 'delta ' + stat )
            pp.suptitle( asm1 + ' vs. ' + asm2 )
            nonzero_deltas=[v for v in deltas if v != 0]
            pp.title( 'no diff for %d of %d (%f %%)' % ( len(deltas)-len(nonzero_deltas),
                                                         len(deltas), 
                                                         asmutil.perc(len(deltas)-len(nonzero_deltas),len(deltas))))

            pp.hist(nonzero_deltas)
            asmutil.ensure_dir_exists( out_dir )
            pp.savefig( os.path.join( out_dir, asmutil.make_alpha_num(stat) + '.png' ) )
            pp.close(fig)

    def import_quast_stats( self, quast_dir, asms ):
        """Import assembly stats computed by QUAST"""
        
        logging.info('quast_dir='+quast_dir)
        for readset in os.listdir(quast_dir):
            logging.info('readset='+readset)
            z = pandas.read_csv( os.path.join( quast_dir, readset, 'report.tsv' ), sep = '\t' )
            for t in z.itertuples(index=False):
                for asm_idx, asm in enumerate(asms):
                    value=t[1+asm_idx]
                    stat=t[0]
                    if asmutil.is_convertible_to(value,float):
                        # logging.info('readset=%s stat=%s asm_idx=%d asm=%s value=%s',
                        #              readset, stat, asm_idx, asm, value)
                        self.put_stat(readset=readset, asm=asm, stat=stat, value=value)

# end: class AsmCmp

def test_AsmCmp():
    logging.basicConfig(level=logging.INFO)
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
                                 readsets = ( 'READS1', 'READS2' ),
                                 stats = [ 'NA50' ], out_dir = os.path.join( statdir, 'plots' ) )
    quast_dir='/idi/sabeti-scratch/ilya/sw/quast-3.2/ebov/assembly-vfat__assembly-mummer3/quast/assembly-vfat__assembly-mummer3'
    readsets=os.listdir(quast_dir)
    stats=list(pandas.read_csv(os.path.join(quast_dir,readsets[0],'report.tsv'),sep='\t',header=0).ix[:,0])
    print('stats=',stats)
    #asmcmp.import_quast_stats( quast_dir=quast_dir, asms='vfat mummer3'.split() )
    logging.info('-------------finding largest changes----------------')
    asmcmp.find_largest_changes( asm1='vfat', asm2='mummer3', readsets=readsets, stats=stats,
                                 out_dir=os.path.join(statdir,'plots2'))

    

if __name__ == '__main__':
    test_AsmCmp()


        
