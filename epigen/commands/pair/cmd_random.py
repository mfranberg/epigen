import click
import random
from math import sqrt

from epigen.commands.command import CommandWithHelp
from epigen.plink import generate

def random_penetrance(H2, p_d, pmin = 0.1, pmax = 0.9):
    return [ min( max( random.normalvariate( p_d, sqrt( H2 * p_d * (1 - p_d) ) ), pmin ), pmax ) for i in range( 9 ) ]

@click.command( 'random', cls = CommandWithHelp, short_help='Samples interaction models.' )
@click.option( '--maf', nargs=2, type=generate.probability, help='Minor allele frequency of the two snps.', default = [0.3, 0.3] )
@click.option( '--sample-size', nargs=2, type=int, help='Number of cases and controls.', default = [2000, 2000] )
@click.option( '--ld', type=generate.probability, help='Strength of LD (ignores second maf).', default = None )
@click.option( '--num-pairs', type=int, help='Number of pairs to generate from each model.', default = 100 )
@click.option( '--num-models', type=int, help='Number of different models to generate.', default = 1 )
@click.option( '--heritability', type=float, help='Approximate heritability of each model.', default = 0.02 )
@click.option( '--base-risk', type=float, help='The base risk of the neutral alleles.', default = 0.5 )
@click.option( '--out', type = click.Path( writable = True ), help='Output plink file.', required = True )
def epigen(maf, sample_size, ld, num_pairs, num_models, heritability, base_risk, out):
    models = [ ( num_pairs, random_penetrance( heritability, base_risk ), 1 ) for i in range( num_models ) ]

    fixed_params = generate.FixedParams( maf, ld, sample_size[ 0 ], sample_size[ 1 ] )
    generate.write_data( fixed_params, models, out )
