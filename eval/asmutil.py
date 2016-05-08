"""Misc assembly-related utils"""

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
