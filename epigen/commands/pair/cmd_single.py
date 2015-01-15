import click
from epigen.commands.command import CommandWithHelp
from epigen.plink import generate

@click.command( 'single', cls = CommandWithHelp, short_help="Generates a plink file by conditioning on the phenotype and generating genotypes, useful for case/control." )
@click.option( '--model', nargs=9, type=generate.probability, help='Space-separated list of floating point numbers that represents the penetrance matrix, specified row-wise from left to right.', required = True )
@click.option( '--maf', nargs=2, type=generate.probability, help='Minor allele frequency of the two snps.', default = [0.3, 0.3] )
@click.option( '--sample-maf/--no-sample-maf', help='The --maf is treated as a range and maf is sampled uniformly in this range.', default = False )
@click.option( '--sample-size', nargs=2, type=int, help='Number of cases and controls.', default = [2000, 2000] )
@click.option( '--npairs', type=int, help='Number of interaction pairs', default = 100 )
@click.option( '--ld', type=generate.probability, help='Strength of LD (ignores second maf).', default = None )
@click.option( '--out', help='Output plink file.', type=click.Path( writable = True ), required = True )
def epigen(model, maf, sample_maf, sample_size, npairs, ld, out):
    fixed_params = generate.FixedParams( maf, ld, sample_size[ 0 ], sample_size[ 1 ], sample_maf )
    models = [ ( npairs, model, False ) ]

    generate.write_data( fixed_params, models, out )
