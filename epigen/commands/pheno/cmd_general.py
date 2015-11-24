import click
from plinkio import plinkfile

from epigen.plink import generate, genmodels, info
from epigen.util import probability
from epigen.commands.command import CommandWithHelp
from epigen.plink.util import find_rows, sample_loci_set, find_beta0, compute_mafs

@click.command( 'general', cls = CommandWithHelp, short_help='Generates a phenotype under the given model and plink file.' )
@click.option( '--model', type=click.Choice( genmodels.get_models( ) ), help='The type of model to use.', required = True )
@click.option( '--mu', nargs=9, type=float, help='Space-separated list of floating point numbers that represents the mean value for each genotype, specified row-wise from left to right.', required = True )
@click.option( '--dispersion', type=float, help='The dispersion parameter (only used in normal for now).', default = 1.0 )
@click.option( '--pair', nargs=2, type=str, help='Name of two SNPs for which the phenotype should be based on (otherwise random).', default = None )
@click.option( "--plink-format/--no-plink-format", help="Use plink format for the phenotype file.", default = False )
@click.option( '--out', type=click.File( "w" ), help='Output phenotype file.', required = True )
@click.argument( 'plink_file', type=click.Path( exists = False ) )
def epigen(model, mu, dispersion, pair, plink_format, out, plink_file):
    input_file = plinkfile.open( plink_file )
    loci = input_file.get_loci( )

    snp_indices = sample_loci_set( loci, 2 )
    if pair:
        snp_indices = [ i for i, x in enumerate( loci ) if x.name in pair ]

    rows = find_rows( input_file, snp_indices )

    mu_map = genmodels.GeneralMuMap( mu )
    pheno_generator = genmodels.get_pheno_generator( model, mu_map, dispersion )
    generate.write_general_phenotype( input_file.get_samples( ), rows, pheno_generator, out, plink_format )
    info.write_info( model, mu, compute_mafs( rows ), dispersion, pheno_generator.sample_size, plink_file + ".info" )
