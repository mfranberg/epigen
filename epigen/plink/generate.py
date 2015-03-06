import random
import os
import sys
from math import exp

from .output import OutputFiles
from .plink_file import PlinkFile
from .genmodels import joint_maf
from epigen.plink import util

from plinkio import plinkfile

##
# Writes the plink data in the location specified by the
# given arguments.
#
# @param m The type of GLM model used to generate data.
# @param fixed_params The simulation parameters.
# @param param_list The list of parameters to generate from.
# @param output_prefix The output prefix, different file endings will be generated.
#
def write_general_data(model, fixed_params, param_list, output_prefix):
    path, ext = os.path.splitext( output_prefix )

    # Number of samples must be known beforehand
    phenotype = model.generate_phenotype( fixed_params )
    output_files = OutputFiles( path, phenotype, model.is_binary( ) )
  
    model_index = 1
    for num_pairs, is_case, params in param_list:
        model.init_cache( fixed_params, params, phenotype )
        for i in range( num_pairs ):
            snp1, snp2 = model.generate_genotype( fixed_params, params, phenotype )
            output_files.write( snp1, snp2, is_case, model_index )
        
        model_index += 1
  
    output_files.close( )
    
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
# Writes the phenotypes for two snps and a set of individuals to a file.
#
# @param sample_list List of plinkio.plinkfile.Sample.
# @param rows List of genotypes for all variants.
# @param beta List of beta coefficients for each variant.
# @param dispersion Dispersion parameter.
# @param link The link function.
# @param output_file The phenotypes will be written to this file.
# @param beta0 The intercept, if not specified will be set to the overall mean.
#
def write_general_phenotype(sample_list, rows, pheno_generator, output_file, plink_format):
    na_string = "NA"
    if plink_format:
        na_string = "-9"

    output_file.write( "FID\tIID\tPheno\n" )
    for i, sample in enumerate( sample_list ):
        variants = [ rows[ j ][ i ] for j in range( len( rows ) ) ]

        pheno = pheno_generator.generate_pheno( variants )
        pheno_str = str( pheno )
        if pheno == None:
            pheno_str = na_string
 
        output_file.write( "{0}\t{1}\t{2}\n".format( sample.fid, sample.iid, pheno_str ) )

##
# Writes the .pair file for a given plink file. This
# can be very time consuming for a large number of variants.
#
def generate_pairs(plink_prefix):
    pf = plinkfile.open( plink_prefix )
    loci = pf.get_loci( )
   
    if len( loci ) > 10000:
        raise ValueError( "Creating pairs for more than 10000 variants is too time consuming." )

    with open( plink_prefix + ".pair", "w" ) as pair_file:
        for i in range( len( loci ) ):
            for j in range( i + 1, len( loci ) ):
                pair_file.write( "{0} {1}\n".format( loci[ i ].name, loci[ j ].name ) )

    return

##
# Generate a set of single variants.
#
def write_single(nvariants, nsamples, output_prefix, maf = None, create_pair = False):
    if create_pair and nvariants > 10000:
        raise ValueError( "Creating pairs for more than 10000 variants is too time consuming." )
    
    pf = PlinkFile( output_prefix, [ 2 ] * nsamples, 0 )

    # These a and b values were taken by fitting a beta distribution to the
    # allele frequency distribution of EUR 1000G.
    generate_maf = lambda: random.betavariate( 0.4679562, 0.4679562 )
    if maf:
        geneate_maf = lambda: maf[ 0 ] + ( maf[ 1 ] - maf[ 0 ] ) * random.random( )
    
    for i in range( nvariants ):
        m = generate_maf( )
        prob = [ (1-m)**2, 2*m*(1-m), m**2 ]
        genotypes = [ util.sample_categorical( prob, [0, 1, 2] ) for j in range( nsamples ) ]
        pf.write( i, genotypes )

    pf.close( )

    if create_pair:
        generate_pairs( output_prefix )

