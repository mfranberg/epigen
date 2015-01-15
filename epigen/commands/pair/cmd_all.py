import click
from functools import partial

from epigen.commands.command import CommandWithHelp
from epigen.plink import generate
from epigen.interaction.generator import InteractionGenerator, mat_or
from epigen.interaction.util import heritability

##
# Given a desired heritability, a fully penetrant disease model
# (penetrance is either 1 or 0), this function tries to find a
# non-fully penetrant disease model that has the desired
# heritability.
#
# @param desired_heritability The desired heritability
# @param base_risk The population risk.
# @param maf Minor allele frequency.
# @param model The fully penetrance disease model.
# 
# @return A non-fully penetrant disease model with as close
#         heritability to the desired as possible.
#
def find_penetrance(desired_heritability, base_risk, maf, model):
    if sum( model ) == 0 or sum( model ) == 9:
        return [ 0.5 ] * 9

    step = 0.01
    disease_p = base_risk + step
    disease_map = dict( { 0 : base_risk, 1 : 1.0 } )
    penetrance = list( map( lambda x: disease_map[ x ], model ) )
    if heritability( penetrance, maf ) < desired_heritability:
        print( "Warning: impossible to find penetrance with desired heritability under (will use maximum): " + str( model ) )

    while disease_p <= 1.0:
        disease_map = dict( { 0 : base_risk, 1 : disease_p } )
        penetrance = list( map( lambda x: disease_map[ x ], model ) )
        if heritability( penetrance, maf ) >= desired_heritability:
            return penetrance

        disease_p += step

    return penetrance

@click.command( 'all', cls = CommandWithHelp, short_help='Generates all possible interaction pairs.' )
@click.option( '--maf', nargs=2, type=generate.probability, help='Minor allele frequency of the two snps.', default = [0.3, 0.3] )
@click.option( '--sample-size', nargs=2, type=int, help='Number of cases and controls.', default = [ 2000, 2000 ] )
@click.option( '--ld', type=generate.probability, help='Strength of LD (ignores second maf).', default = None )
@click.option( '--num-pairs', type=int, help='Number of pairs to generate from each model.', default = 100 )
@click.option( '--heritability', type=float, help='Approximate heritability of each model.', default = 0.02 )
@click.option( '--base-risk', type=float, help='The base risk of the neutral alleles.', default = 0.5 )
@click.option( '--out', type = click.Path( writable = True ), help='Output .tped file.', required = True )
def epigen(maf, sample_size, ld, num_pairs, heritability, base_risk, out):
    models = [ ]
    generator = InteractionGenerator( mat_or )
    interactions, nulls = generator.generate( )

    partial_find_penetrance = partial( find_penetrance, heritability, base_risk, maf )
    interaction_penetrances = list( map( partial_find_penetrance, interactions ) )

    models = [ ( num_pairs, p, 1 ) for p in interaction_penetrances ]

    fixed_params = generate.FixedParams( maf, ld, sample_size[ 0 ], sample_size[ 1 ] )
    generate.write_data( fixed_params, models, out )
