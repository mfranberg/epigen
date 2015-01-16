import random

##
# Given a list of variant names this function finds
# the corresponding genotypes and returns them.
#
def find_rows(plink_file, loci):
    loci_set = set( loci )
    rows = [ ]
    for i, row in enumerate( plink_file ):
        if i in loci_set:
            rows.append( list( row ) )

    return rows

##
# Returns the index of the given variant names.
#
# @param loci List of loci.
# @param names A list of variant names.
#
# @return The indices of the given locus.
#
def find_index(loci, names):
    loci_name_set = set( names )
    return [ i for i, x in enumerate( loci ) if x in loci_name_set ]

##
# Randomly selects n loci from the loci list.
#
# @param loci List of loci.
#
# @return The indices of the selected loci.
#
def sample_loci_set(loci, n):
    loci_index = list( range( len( loci ) ) )
    random.shuffle( loci_index )

    return loci_index[:n]
