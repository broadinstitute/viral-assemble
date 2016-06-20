"""Misc assembly-related utils"""

from __future__ import division
import os, re

def dump_file( fname, value ):
    """store string in file"""
    with open( fname, 'w' )  as out: 
        out.write( str( value ) )

def slurp_file( fname ):
    """Read entire file into one string"""
    if not os.path.isfile( fname ): raise IOError( 'File not found: %s' % fname )
    with open( fname ) as f:
        return f.read()

        
def ensure_dir_exists( d ):
    """Create the given dir (and any parent dirs) if it does not already exist"""
    if not os.path.isdir( d ): os.makedirs( d )

def make_alpha_num( s ):
    """Return a version of the argument string, in which all non-alphanumeric chars have been replaced
    by underscores.
    """
    return re.sub( '\W+', '_', s )
    
def is_convertible_to( x, aType ):
    """Test if a value is convertible to a given type"""
    try:
        dummy = aType(x)
        return True
    except ValueError: return False

def perc(num,den):
    """Compute an integer percentage"""
    assert float(num)>=0 and float(den)>0
    return int(100.0*(float(num)/float(den)))

def dict_get( d, k, dflt ):
    """Get a value from the dictionary; if not such value, return the specified default"""
    return d.get( k, dflt )


def tabjoin( *args, **kwargs ):
    """Join args by tab"""
    sep = dict_get( kwargs, 'sep', '\t' )
    return sep.join( map( str, args ) )

def tabwriten( f, *args, **kwargs ):
    """Write a tab-separated line to a file"""
    f.write( tabjoin( *args, **kwargs ) )

def tabwrite( f, *args, **kwargs ):
    """Write a tab-separated line to a file, followed by a newline"""
    tabwriten( f, *args, **kwargs )
    f.write( '\n' )
    
