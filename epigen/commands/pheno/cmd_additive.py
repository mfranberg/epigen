import click
from plinkio import plinkfile
import random

from epigen.plink import generate, genmodels, info
from epigen.plink.util import find_rows, sample_loci_set, find_beta0, generate_beta, compute_mafs
from epigen.commands.command import CommandWithHelp

@click.command( 'additive', cls = CommandWithHelp, short_help='Generates binary phenotypes for given plink data.' )
@click.argument( 'plink_file', type=click.Path( ) )
@click.option( '--beta0', type=float, help='Sets the intercept, by default it is chosen to get 50/50 cases and controls.', default = None )
@click.option( '--beta', nargs=2, type=float, help='The mean and variance of the beta variables (taken from a normal).', required = True )
@click.option( '--num-loci', type=int, help='The number of loci that is involved in the phenotype.', default = 10 )
@click.option( '--model', type=click.Choice( genmodels.get_models( ) ), help="The model to use.", required = True )
@click.option( '--link', type=click.Choice( genmodels.get_links( ).keys( ) ), help="The link function to use.", default = "default" )
@click.option( '--dispersion', type=float, help="The dispersion parameter to use.", default=1.0 )
@click.option( '--out', type = click.File( 'w' ), help='Output phenotype file.', required = True )
def epigen(plink_file, beta0, beta, num_loci, model, link, dispersion, out):
    input_file = plinkfile.open( plink_file ) 
    loci = input_file.get_loci( )
    snp_indices = sample_loci_set( loci, num_loci )
    gen_beta = generate_beta( num_loci, beta[ 0 ], beta[ 1 ] )
    rows = find_rows( input_file, snp_indices )
    
    with open( plink_file + ".av", "w" ) as av_file:
        for index, b in zip( snp_indices, gen_beta ):
            av_file.write( loci[ index ].name + " " + str( b ) + "\n" )

    if not beta0:
        beta0 = find_beta0( rows, gen_beta )

    mu_map = genmodels.AdditiveMuMap( beta0, gen_beta, genmodels.get_link( model, link ) )
    pheno_generator = genmodels.get_pheno_generator( model, mu_map, dispersion )
    generate.write_general_phenotype( input_file.get_samples( ), rows, pheno_generator, out, False )
    info.write_info( model, mu_map, compute_mafs( rows ), dispersion, pheno_generator.sample_size, plink_file + ".info" )
