import click
from epigen.commands.command import CommandWithHelp
from epigen.plink import generate, genmodels
from epigen.util import probability

@click.command( 'glm', cls = CommandWithHelp, short_help="Generates a plink file by conditioning on the phenotype and generating genotypes, useful for case/control." )
@click.option( '--model', type=click.Choice( genmodels.get_models( ) ), help = "The distribution of the phenotype around the mean.", required = True )
@click.option( '--link', type=click.Choice( genmodels.get_links( ).keys( ) ), help='The link function to use', required = True )
@click.option( '--beta', nargs=9, type=float, help='Space-separated list of regression coefficients a, b1, b2, g1, g2, d11, d12, d21 and d22.', required = True )
@click.option( '--dispersion', nargs=1, type=float, help='Dispersion parameter (only used for normal atm).', default = 1.0 )
@click.option( '--maf', nargs=2, type=probability.probability, help='Minor allele frequency of the two snps.', default = [0.3, 0.3] )
@click.option( '--sample-maf/--no-sample-maf', help='The --maf is treated as a range and maf is sampled uniformly in this range.', default = False )
@click.option( '--sample-size', nargs=2, type=int, help='Number of samples (for binomial cases and controls, only first will be considered otherwise).', default = [2000, 2000] )
@click.option( '--npairs', type=int, help='Number of interaction pairs', default = 100 )
@click.option( '--ld', type=probability.probability, help='Strength of LD (ignores second maf).', default = None )
@click.option( '--out', help='Output plink file.', type=click.Path( writable = True ), required = True )
def epigen(model, link, beta, dispersion, maf, sample_maf, sample_size, npairs, ld, out):
    lf = genmodels.get_link( link )

    mu = genmodels.get_mean_values( beta, lf )
    model_def, params = genmodels.get_model_and_params( model, mu, dispersion, maf )

    fixed_params = genmodels.FixedParams( maf, ld, sample_size, sample_maf )
    model_list = [ ( npairs, 1, params ) ]
    generate.write_general_data( model_def, fixed_params, model_list, out )
