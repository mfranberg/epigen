import numpy as np
import argparse
import random
import os
from math import sqrt, exp

from plinkio import plinkfile

from .output import OutputFiles

##
# The parameters that does not changed between models.
# 
class FixedParams:
    def __init__(self, maf, ld, ncases, ncontrols, sample_maf = False):
        self.maf = maf
        self.ld = ld
        self.ncases = ncases
        self.ncontrols = ncontrols
        self.sample_maf = sample_maf

    def get_maf(self):
        if not self.sample_maf:
            return self.maf
        else:
            start = self.maf[ 0 ]
            end = self.maf[ 1 ] - start

            return start + end * np.random.random( 2 )

##
# Converts the given value to a floating point value
# if it represents a probability, otherwise it raises
# a ValueError.
#
# @param prob_string String representing a probability.
#
# @return String as a floating point probability.
#
def probability(prob_string):
    prob_value = float( prob_string )
    if 0.0 <= prob_value <= 1.0:
        return prob_value
    else:
        raise ValueError( "Probability not between 0 and 1." )

    return prob_value

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
# Computes the probability of a pair of snps x1, x2
# given the disease status.
#
# @param penetrance Penetrance of each genotype.
# @param maf Desired minor allele frequncy.
# @param ld The probability Pr[x2 = i| x1 = i].
#
# @return A tuple containing first the probabilities of 
#         x1, x2 | D = 1, and second the probabilityies
#         of x1, x2 | D = 2 represented as vectors.
#
def joint_snp(penetrance, maf, ld):
    joint_hw = np.array( joint_maf( maf, ld ) )

    penetrance = np.array( penetrance )

    joint_numerator = joint_hw * penetrance
    denom = sum( joint_numerator )
    sick_prob = ( joint_numerator / denom ).tolist( )

    joint_numerator = joint_hw * ( 1 - penetrance )
    denom = sum( joint_numerator )
    healthy_prob = ( joint_numerator / denom ).tolist( )

    return ( sick_prob, healthy_prob )

##
# Samples the given number of samples from
# the joint probability model.
#
# @param joint_prob Joint probability model represented as 
#                   a vector.
# @param num_samples Number of samples to take.
#
# @return Iterator over all samples.
#
def sample(joint_prob, num_samples):
    samples = np.random.multinomial( num_samples, joint_prob, 1 )[ 0 ]
    for i in range( 3 ):
        for j in range( 3 ):
            for k in range( samples[ 3 * i + j ] ):
                yield (i, j)

##
# Generates two lists of genotypes for two SNPs under
# the given joint probability model.
#
# @param joint_prob The probability model.
# @param num_samples Number of samples.
#
# @return A tuple containing first the genotypes of the
#         first snp, and second the genotypes of the second
#         snp.
#
def generate_genotypes(joint_prob, num_samples):
    snp1_list = [ 0 ] * num_samples
    snp2_list = [ 0 ] * num_samples
    i = 0
    snps = [ (snp1, snp2) for snp1, snp2 in sample( joint_prob, num_samples ) ]
    random.shuffle( snps )

    for snp1, snp2 in snps:
        snp1_list[ i ] = snp1
        snp2_list[ i ] = snp2

        i += 1

    return ( snp1_list, snp2_list )

##
# Writes the plink data to the given .tped, .case and .pair files.
#
# @param num_pairs The number of pairs to generate.
# @param penetrance The model penetrance.
# @param is_case Determines whether this is an interaction or not.
# @param fixed_params The simulation parameters.
# @param model_index Index of the model.
# @param output_files The output files.
#
def write_genotypes(num_pairs, penetrance, is_case, fixed_params, model_index, output_files):
    for i in range( num_pairs ):
        maf = fixed_params.get_maf( )
        
        sick_prob, healthy_prob = joint_snp( penetrance, maf, fixed_params.ld )
        snp1_sick, snp2_sick = generate_genotypes( sick_prob, fixed_params.ncases )
        snp1_healthy, snp2_healthy = generate_genotypes( healthy_prob, fixed_params.ncontrols )
   
        output_files.write( snp1_sick + snp1_healthy, snp2_sick + snp2_healthy, is_case, model_index )

##
# Writes the plink data in the location specified by the
# given arguments.
#
# @param fixed_params The simulation parameters.
# @param models The list of models to generate from.
# @param output_prefix The output prefix, different file endings will be generated.
#
def write_data(fixed_params, models, output_prefix):
    path, ext = os.path.splitext( output_prefix )

    # Number of samples must be known beforehand
    output_files = OutputFiles( path, fixed_params.ncases, fixed_params.ncontrols )
  
    model_index = 1
    for num_pairs, params, is_case in models:
        write_genotypes( num_pairs, params, is_case, fixed_params, model_index, output_files )
        model_index += 1
  
    output_files.close( )
    
##
# Generates two distinct loci.
#
# @param loci List of locus.
#
# @return A list of two loci.
#
def sample_two_loci(loci):
    loci_index = list( range( len( loci ) ) )
    random.shuffle( loci_index )

    return loci_index[:2]

##
# Generates a phenotype for the given snp pair and penetrance
# matrix. 
#
# @param snp1 Genotype of first snp
# @param snp2 Genotype of second snp
# @param model Penetrance matrix as a length 9 list
#
# @return 1 or 0 representing case control if no snp was missing, None
#         otherwise.
#
def generate_phenotype(snp1, snp2, model):
    if snp1 != 3 and snp2 != 3:
        return int( random.random( ) <= model[ 3 * snp1 + snp2 ] )
    else:
        return None

##
# Generates a normal phenotype for the given snp pair and mean
# matrix. 
#
# @param snp1 Genotype of first snp
# @param snp2 Genotype of second snp
# @param model Mean matrix as a length 9 list
# @param sd Common standard deviation
#
# @return Normal variable if not missing, None otherwise.
#
def generate_phenotype_cont(snp1, snp2, model, sd):
    if snp1 != 3 and snp2 != 3:
        return random.normalvariate( model[ 3 * snp1 + snp2 ], sd )
    else:
        return None

##
# Generates a phenotype for the given snp pair and penetrance
# matrix. 
#
# @param snp1 Genotype of first snp
# @param snp2 Genotype of second snp
# @param model Penetrance matrix as a length 9 list
#
# @return 1 or 0 representing case control if no snp was missing, None
#         otherwise.
#
def generate_phenotype_additive(variants, beta0, beta):
    if 3 in variants:
        return None

    score = beta0 + sum( [ b * v for b, v in zip( beta, variants ) ] )
    p = 1/(1 + exp(-score))
    
    return int( random.random( ) <= p )

##
# Generates a phenotype for the given snp pair and penetrance
# matrix. 
#
# @param snp Genotype of variant.
# @param env Environmental factor.
# @param model Mean for each genotype and environmental level.
# @param std Standard deviation of phenotype.
#
# @return A continuous phenotype, or "NA" if genotype was missing.
#
def generate_env_phenotype(snp, env, model, std):
    if snp != 3:
        return random.normalvariate( model[ env * 3 + snp ], std )
    else:
        return "NA"

##
# Generates an environment variable
# from the level frequencies.
#
# @param level frequencies that sum to one.
#
# @return The generated level.
#
def generate_environment(env):
    u = random.random( )
    cumsum = 0.0
    for i, e in enumerate( env ):
        cumsum += e
        if u <= cumsum:
            return i

    raise Exception( "Environmental frequencies does not sum to 1." )
