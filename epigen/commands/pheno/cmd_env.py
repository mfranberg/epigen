import click
from plinkio import plinkfile

from epigen.plink import generate, util
from epigen.commands.command import CommandWithHelp

##
# Simple class for storing the parameters
# of the simulation.
#
class EnvParams:
    ##
    # @param env The environment factor.
    # @param model The mean level for each combination.
    # @param std The standard deviation.
    #
    def __init__(self, env, model, std):
        self.env = env
        self.model = model
        self.std = std

##
# Converts a comma-separated string of floats to
# a list of floats.
#
# @param model_str A comma-separated string of floats.
#
# @return A list of floats.
#
def model_type(model_str):
    return [ float( x ) for x in model_str.strip( ).split( "," ) ]

##
# Writes the phenotypes for two snps and a set of individuals to a file.
#
# @param sample_list List of plinkio.plinkfile.Sample.
# @param snp_row List of genotypes for the first snp.
# @param params Parameters for generation
# @param env_file Output file for environmental factors.
# @param pheno_file Output file for phenotype.
#
def write_phenotypes(sample_list, snp_row, params, env_file, pheno_file):
    env_file.write( "FID\tIID\tEnv\n" )
    pheno_file.write( "FID\tIID\tPheno\n" )

    for sample, snp in zip( sample_list, snp_row ):
        env = generate.generate_environment( params.env )
        pheno = generate.generate_env_phenotype( snp, env, params.model, params.std )

        env_file.write( "{0}\t{1}\t{2}\n".format( sample.fid, sample.iid, str( env ) ) )
        pheno_file.write( "{0}\t{1}\t{2}\n".format( sample.fid, sample.iid, str( pheno ) ) )

@click.command( 'env', cls = CommandWithHelp, short_help='Generates a phenotype and environment variable for a gene-environment model.' )
@click.option( '--model', type=model_type, help='Comma-separated list of mean for each cell, columns correspond to the genotype.', required = True )
@click.option( '--std', help='Common standard deviation.', default = 1.0 )
@click.option( '--env', type=model_type, help='Frequency of each environmental level, must sum to 1, default is equal frequency.', default = None )
@click.option( '--variant', type=str, help='Name of the variant the phenotype should be based on, default is random.', default = None )
@click.option( '--pheno-file', help='Output phenotype file.', type = click.File( "w" ), required = True )
@click.option( '--env-file', help='Output environment file.', type = click.File( "w" ), required = True )
@click.argument( 'plink_file', type=click.Path( writable=True ) )
def epigen(model, std, env, variant, pheno_file, env_file, plink_file):
    input_file = plinkfile.open( plink_file )
    loci = input_file.get_loci( )

    snp_index = util.sample_loci_set( input_file.get_loci( ), 1 )[ 0 ]
    if variant:
        snp_index = util.find_index( loci, variant )

    snp_row = util.find_rows( input_file, [ snp_index ] )[ 0 ]

    env_freq = env
    if not env_freq:
        num_env = len( args.model ) / 3
        env_freq = [ 1.0 / num_env ] * num_env

    params = EnvParams( env_freq, model, std )
    write_phenotypes( input_file.get_samples( ), snp_row, params, env_file, pheno_file )
