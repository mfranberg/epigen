import random

##
# Given a list of loci nmaes this function finds
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
