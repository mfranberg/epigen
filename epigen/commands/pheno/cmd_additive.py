import click
from plinkio import plinkfile
import random

from epigen.plink import generate, genmodels
from epigen.plink.util import find_rows, sample_loci_set, find_beta0
from epigen.commands.command import CommandWithHelp

@click.command( 'additive', cls = CommandWithHelp, short_help='Generates binary phenotypes for given plink data.' )
@click.argument( 'plink_file', type=click.Path( ) )
@click.option( '--beta0', type=float, help='Sets the intercept, by default it is chosen to get 50/50 cases and controls.', default = None )
@click.option( '--beta', nargs=2, help='The mean and variance of the beta variables (taken from a normal).', required = True )
@click.option( '--num-loci', type=int, help='The number of loci that is involved in the phenotype.', default = 10 )
@click.option( '--model', type=click.Choice( genmodels.get_models( ) ), help="The model to use.", required = True )
@click.option( '--link', type=click.Choice( genmodels.get_links( ).keys( ) ), help="The link function to use.", required = True )
@click.option( '--dispersion', type=float, help="The dispersion parameter to use.", default=1.0 )
@click.option( '--out', type = click.Path(writable=True), help='Output phenotype file.', required = True )
def epigen(plink_file, beta0, beta, num_loci, model, link, dispersion, out):
    input_file = plinkfile.open( plink_file ) 
    loci = input_file.get_loci( )
    snp_indices = sample_loci_set( loci, num_loci )
    gen_beta = generate_beta( num_loci, beta[ 0 ], beta[ 1 ] )
    rows = find_rows( input_file, snp_indices )

    if not beta0:
        beta0 = find_beta0( rows, gen_beta )

    mu_map = genmodels.AdditiveMuMap( beta0, gen_beta, genmodels.get_link( link ) )
    pheno_generator = genmodels.get_pheno_generator( model, mu_map, dispersion )
    generate.write_general_phenotype( input_file.get_samples( ), rows, pheno_generator, out, plink_format )
