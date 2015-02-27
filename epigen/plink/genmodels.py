import random
from math import sqrt

import numpy as np
from scipy.stats import norm

##
# The parameters that does not changed between models.
# 
class FixedParams:
    def __init__(self, maf, ld, sample_size, sample_maf = False):
        self.maf = maf
        self.ld = ld
        self.sample_size = sample_size
        self.sample_maf = sample_maf

    def get_maf(self):
        if not self.sample_maf:
            return joint_maf( self.maf, self.ld )
        else:
            start = self.maf[ 0 ]
            end = self.maf[ 1 ] - start

            return joint_maf( start + end * np.random.random( 2 ), self.ld )

    def num_samples(self):
        return self.sample_size

    def num_cases(self):
        return self.sample_size[ 0 ]

    def num_controls(self):
        return self.sample_size[ 1 ]

##
# This class represents a binary model that is used to first generate
# a phenotype, and then generate the genotypes.
#
class BinaryModel:
    def generate_phenotype(self, fixed_params):
        return [ 1 ] * fixed_params.num_cases( ) + [ 0 ] * fixed_params.num_controls( )

    def joint_prob(self, maf, penetrance, phenotype):
        geno_prob = [ p * f for p, f in zip( penetrance, maf ) ]
        geno_denom = sum( geno_prob )

        return [ g / geno_denom for g in geno_prob ]

    def generate_genotype(self, fixed_params, params, phenotype):
        maf = fixed_params.get_maf( )
        prob = [ ]
        prob.append( self.joint_prob( maf, params.penetrance, 1 ) )
        prob.append( self.joint_prob( maf, params.penetrance, 0 ) )

        snp1_list = list( )
        snp2_list = list( )
        for pheno in phenotype:
            snp1, snp2 =  sample_categorical( prob[ pheno ] )
            snp1_list.append( snp1 )
            snp2_list.append( snp2 )

        return snp1_list, snp2_list

    def is_binary(self):
        return True

class BinaryParams:
    def __init__(self, penetrance):
        self.penetrance = penetrance
    
##
# This class represents a continuous model that is used to first generate
# a phenotype, and then generate the genotypes.
#
class ContModel:
    def __init__(self, mu, std, maf):
        self.mu = mu
        self.std = std
        self.maf = joint_maf( maf, None )

    def generate_phenotype(self, fixed_params):
        mu = sum( ( m * f for m, f in zip( self.mu, self.maf ) ) )
        mu2 = sum( ( m**2 * f for m, f in zip( self.mu, self.maf ) ) )
        std2 = sum( ( s**2 * f for s, f in zip( self.std, self.maf ) ) )

        total_std = sqrt( mu2 - mu**2 + std2 )

        return [ random.normalvariate( mu, total_std ) for i in range( fixed_params.num_samples( ) ) ]

    def joint_prob(self, mu, std, maf, pheno):
        geno_prob = [ norm.pdf( pheno, m, s ) * f for m, s, f in zip( mu, std, maf ) ]
        geno_denom = sum( geno_prob )

        return [ g / geno_denom for g in geno_prob ]

    def generate_genotype(self, fixed_params, params, phenotype):
        snp1_list = list( )
        snp2_list = list( )
        maf = fixed_params.get_maf( )
        for pheno in phenotype:
            prob_geno = self.joint_prob( params.mu, params.std, maf, pheno )
            snp1, snp2 = sample_categorical( prob_geno )
            
            snp1_list.append( snp1 )
            snp2_list.append( snp2 )

        return snp1_list, snp2_list
    
    def is_binary(self):
        return False

class ContParams:
    def __init__(self, mu, std):
        self.mu = mu
        self.std = std

##
# Sample genotyeps from a categorical distribution.
#
# @param prob The probability for each genotype, specified as a vector by row.
#
# @return The sampled genotype.
#
def sample_categorical(prob):
    r = random.random( )
    cum = 0.0
    for c, p in zip( [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)], prob ):
        cum += p
        if r <= cum:
            return c

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
