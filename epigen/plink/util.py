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

##
# Sample genotyeps from a categorical distribution.
#
# @param prob The probability for each genotype, specified as a vector by row.
#
# @return The sampled genotype.
#
def sample_categorical(prob, cat=[(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]):
    r = random.random( )
    cum = 0.0
    for c, p in zip( cat, prob ):
        cum += p
        if r <= cum:
            return c

def fast_sample_setup(probs):
    K = len( probs )
    q = [ 0.0 ] * K
    J = [ 0 ] * K
 
    # Sort the data into the outcomes with probabilities
    # that are larger and smaller than 1/K.
    smaller = [ ]
    larger  = [ ]
    for kk, prob in enumerate( probs ):
        q[ kk ] = K * prob
        if q[ kk ] < 1.0:
            smaller.append( kk )
        else:
            larger.append( kk )
 
    # Loop though and create little binary mixtures that
    # appropriately allocate the larger outcomes over the
    # overall uniform mixture.
    while len( smaller ) > 0 and len( larger ) > 0:
        small = smaller.pop( )
        large = larger.pop( )
 
        J[ small ] = large
        q[ large ] = q[ large ] - ( 1.0 - q[ small ] )
 
        if q[ large ] < 1.0:
            smaller.append( large )
        else:
            larger.append( large )
 
    return J, q
 
def fast_sample(J, q, cat=[(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)] ):
    K  = len( J )
 
    # Draw from the overall uniform mixture.
    kk = random.randint( 0, K - 1 )
 
    # Draw from the binary mixture, either keeping the
    # small one, or choosing the associated larger one.
    if random.random( ) < q[ kk ]:
        return cat[ kk ]
    else:
        return cat[ J[ kk ] ]

##
# Computes the joint Hardy-Weinberg model represented
# as a vector.
#
# @param maf Desired minor allele frequncy.
# @param ld The probability Pr[x2 = i| x1 = i].
#
# @return Joint probability model for the Hardy-Weinberg
#         model.
#
def joint_maf(maf, ld):
    maf1 = maf[ 0 ]
    maf2 = maf[ 1 ]

    p = [ ( 1 - maf[ 0 ] )**2, 2 * maf[ 0 ] * ( 1 - maf[ 0 ] ), ( maf[ 0 ] )**2 ]
    q = [ ( 1 - maf[ 1 ] )**2, 2 * maf[ 1 ] * ( 1 - maf[ 1 ] ), ( maf[ 1 ] )**2 ]

    if not ld:
        return   [ p[ 0 ] * q[ 0 ], p[ 0 ] * q[ 1 ], p[ 0 ] * q[ 2 ],
                   p[ 1 ] * q[ 0 ], p[ 1 ] * q[ 1 ], p[ 1 ] * q[ 2 ],
                   p[ 2 ] * q[ 0 ], p[ 2 ] * q[ 1 ], p[ 2 ] * q[ 2 ] ]
    else:
        return   [ p[ 0 ] * ld, p[ 0 ] * ( 1 - ld ) / 2.0, p[ 0 ] * ( 1 - ld ) / 2.0,
                   p[ 1 ] * ( 1 - ld ) / 2.0, p[ 1 ] * ld, p[ 1 ] * ( 1 - ld ) / 2.0,
                   p[ 2 ] * ( 1 - ld ) / 2.0, p[ 2 ] * ( 1 - ld ) / 2.0, p[ 2 ] * ld ]

    return joint_prob

##
# Computes the second allele frequency for the
# given genotypes.
#
# @param row Genotypes 0, 1, 2. Here 3 is missing.
#
# @return The second allele frequency.
#
def compute_maf(row):
    no_missing = filter( lambda x: x != 3, row )
    p = sum( no_missing ) / ( 2.0 * len( no_missing ) )

    return p

##
# Determines a beta0 that gives the probability of being
# a case 0.5.
#
# @param rows List of genotypes for each locus.
# @param beta The beta used for each locus.
#
# @return A beta0 that makes the probability of being a case 0.5.
#
def find_beta0(rows, beta):
    maf = [ compute_maf( r ) for r in rows ]
    beta0 = -sum( b * 2 * m for b, m in zip( beta, maf ) )

    return beta0

##
# Generates effect sizes according to the effect size
# distribution.
#
# @param n The number of beta
# @param mean The mean of the beta distribution
# @param sd The standard deviation of the beta distibution
#
# @return list of effect sizes
#
def generate_beta(n, mean, sd):
    return [ random.normalvariate( mean, sd ) for i in range( n ) ]


