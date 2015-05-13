import click
from epigen.commands.command import CommandWithHelp
from epigen.plink import generate, genmodels, info
from epigen.util import probability

@click.command( 'general', cls = CommandWithHelp, short_help="Generates a plink file by conditioning on the phenotype and generating genotypes, useful for case/control." )
@click.option( '--model', type=click.Choice( genmodels.get_models( ) ), help='The type of model to use.', required = True )
@click.option( '--mu', nargs=9, type=float, help='Space-separated list of floating point numbers that represents the mean value for each genotype, specified row-wise from left to right.', required = True )
@click.option( '--dispersion', type=float, help='The dispersion parameter (only used in normal for now).', default = 1.0 )
@click.option( '--maf', nargs=2, type=probability.probability, help='Minor allele frequency of the two snps.', default = [0.4, 0.4] )
@click.option( '--sample-maf/--no-sample-maf', help='The --maf is treated as a range and maf is sampled uniformly in this range.', default = False )
@click.option( '--sample-size', nargs=2, type=int, help='Number of samples (if only one group only first argument will be used).', default = [2000, 2000] )
@click.option( '--npairs', type=int, help='Number of interaction pairs', default = 100 )
@click.option( '--ld', type=probability.probability, help='Strength of LD (ignores second maf).', default = None )
@click.option( '--out', help='Output plink file.', type=click.Path( writable = True ), required = True )
def epigen(model, mu, dispersion, maf, sample_maf, sample_size, npairs, ld, out):
    fixed_params = genmodels.FixedParams( maf, ld, sample_size, sample_maf )

    model_def, model_params = genmodels.get_model_and_params( model, mu, dispersion, maf )
    params = [ ( npairs, 1, model_params ) ]
    info.write_info( model, mu, maf, dispersion, sample_size, out + ".info" )
    generate.write_general_data( model_def, fixed_params, params, out )
