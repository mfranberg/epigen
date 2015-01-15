import click

from epigen.plink import generate
from epigen.commands.command import CommandWithHelp

##
# Parses a set of models in a given model file.
#
# @param model_file A file where each row contains: number of pairs,
#                   boolean, and 9 penetrance values.
#
# @return An iterator over the models in the file.
#
def parse_models(model_file):
    for line in model_file:
        columns = line.strip( ).split( )

        if len( columns ) != 11:
            print( "Could not parse line in model file." )
            exit( 1 )

        try:
            num_pairs = int( columns[ 0 ] )
            is_case = bool( int( columns[ 1 ], 2 ) )
            params = list( map( float, columns[ 2: ] ) )

            yield (num_pairs, params, is_case)

        except ValueError:
            print( "Could not parse value in model file." )
            exit( 1 )

@click.command( "mixed", cls = CommandWithHelp, short_help="Generates a plink file that contains variant pairs that are generated from different interaction models." )
@click.option( '--model-file', type=click.File( "r" ), help='File where each row contains the number of pairs and 9 penetrance values.', required = True )
@click.option( '--maf', nargs=2, type=generate.probability, help='Minor allele frequency of the two snps.', default = [0.3, 0.3] )
@click.option( '--sample-size', nargs=2, type=int, help='Number of cases and controls', default = [2000, 2000] )
@click.option( '--ld', type=generate.probability, help='Strength of LD (ignores second maf).', default = None )
@click.option( '--out', type = click.Path( writable = True ), help='Output plink file.', required = True )
def epigen(model_file, maf, sample_size, ld, out):
    fixed_params = generate.FixedParams( maf, ld, sample_size[ 0 ], sample_size[ 1 ] )
    models = parse_models( model_file )
    generate.write_data( fixed_params, models, out )
