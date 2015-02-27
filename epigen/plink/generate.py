import os
from math import exp

import numpy as np

from .output import OutputFiles
from .plink_file import PlinkFile
from .genmodels import joint_maf

##
# Writes the plink data in the location specified by the
# given arguments.
#
# @param m The type of GLM model used to generate data.
# @param fixed_params The simulation parameters.
# @param models The list of models to generate from.
# @param output_prefix The output prefix, different file endings will be generated.
#
def write_data(m, fixed_params, models, output_prefix):
    path, ext = os.path.splitext( output_prefix )

    # Number of samples must be known beforehand
    phenotype = m.generate_phenotype( fixed_params )
    output_files = OutputFiles( path, phenotype, m.is_binary( ) )
  
    model_index = 1
    for num_pairs, is_case, params in models:
        for i in range( num_pairs ):
            snp1, snp2 = m.generate_genotype( fixed_params, params, phenotype )
            output_files.write( snp1, snp2, is_case, model_index )
        
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


##
# Generate a set of single variants.
#
def write_single(nvariants, nsamples, output_prefix, maf = None):
    pf = PlinkFile( output_prefix, nsamples, 0 )

    generate_maf = lambda: np.random.beta( 0.8, 0.8 )
    if maf:
        geneate_maf = lambda: maf[ 0 ] + ( maf[ 1 ] - maf[ 0 ] ) * random.random( )
    
    for i in range( nvariants ):
        m = generate_maf( )
        genotypes = np.random.binomial( 2, m, nsamples )
        pf.write( i, genotypes )

    pf.close( )
