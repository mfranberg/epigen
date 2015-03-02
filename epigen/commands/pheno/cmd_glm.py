import click
from plinkio import plinkfile

from epigen.plink import generate, genmodels
from epigen.util import probability
from epigen.commands.command import CommandWithHelp
from epigen.plink.util import find_rows, sample_loci_set, find_beta0

@click.command( 'glm', cls = CommandWithHelp, short_help='Generates a phenotype under the given GLM model and plink file.' )
@click.option( '--model', type=click.Choice( genmodels.get_models( ) ), help='The type of model to use.', required = True )
@click.option( '--link', type=click.Choice( genmodels.get_links( ).keys( ) ), help='The link function to use', default = "default" )
@click.option( '--beta', nargs=9, type=float, help='Space-separated list of regression coefficients a, b1, b2, g1, g2, d11, d12, d21 and d22.', required = True )
@click.option( '--dispersion', type=float, help='The dispersion parameter (only used in normal for now).', default = 1.0 )
@click.option( '--pair', nargs=2, type=str, help='Name of two SNPs for which the phenotype should be based on (otherwise random).', default = None )
@click.option( "--plink-format/--no-plink-format", help="Use plink format for the phenotype file.", default = False )
@click.option( '--out', type=click.File( "w" ), help='Output phenotype file.', required = True )
@click.argument( 'plink_file', type=click.Path( exists = False ) )
def epigen(model, link, beta, dispersion, pair, plink_format, out, plink_file):
    input_file = plinkfile.open( plink_file )
    loci = input_file.get_loci( )

    snp_indices = sample_loci_set( loci, 2 )
    if pair:
        snp_indices = [ i for i, x in enumerate( loci ) if x in pair ]

    rows = find_rows( input_file, snp_indices )
    
    lf = genmodels.get_link( model, link )
    mu = genmodels.get_mean_values( beta, lf )
    
    mu_map = genmodels.GeneralMuMap( mu )
    pheno_generator = genmodels.get_pheno_generator( model, mu_map, dispersion )
    generate.write_general_phenotype( input_file.get_samples( ), rows, pheno_generator, out, plink_format )
