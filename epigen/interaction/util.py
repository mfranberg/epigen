##
# Computes the heritability of a given penetrance and minor
# allele frequency.
#
# @param penetrance 9-element vector of disease penetrance.
# @param maf 2-element vector of minor allele frequencies.
#
# @return The narrow sense heritability.
#
def heritability(penetrance, maf):
    maf1 = maf[ 0 ]
    maf2 = maf[ 1 ]

    p = [ ( 1 - maf[ 0 ] )**2, 2 * maf[ 0 ] * ( 1 - maf[ 0 ] ), ( maf[ 0 ] )**2 ]
    q = [ ( 1 - maf[ 1 ] )**2, 2 * maf[ 1 ] * ( 1 - maf[ 1 ] ), ( maf[ 1 ] )**2 ]

    joint_maf =  [ p[ 0 ] * q[ 0 ], p[ 0 ] * q[ 1 ], p[ 0 ] * q[ 2 ],
                   p[ 1 ] * q[ 0 ], p[ 1 ] * q[ 1 ], p[ 1 ] * q[ 2 ],
                   p[ 2 ] * q[ 0 ], p[ 2 ] * q[ 1 ], p[ 2 ] * q[ 2 ] ]

    pop_p = sum( p * m for p, m in zip( penetrance, joint_maf ) )

    h = 0.0
    for i in range( 3 ):
        for j in range( 3 ):
            cell = 3 * i + j
            h += ( penetrance[ cell ] - pop_p )**2 * joint_maf[ cell ]

    return h / ( pop_p * ( 1 - pop_p ) )
