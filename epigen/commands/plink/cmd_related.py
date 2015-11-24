import click
from epigen.commands.command import CommandWithHelp
from epigen.plink import generate
from epigen.util import probability

@click.command( 'data', cls = CommandWithHelp, short_help="Generates a plink file without a phenotype with related individuals." )
@click.option( '--maf', nargs=2, type=probability.probability, help='If set MAF is generated uniformly between these two values (default uses beta distribution estimated from 1000G).', default = None )
@click.option( '--nsamples', type=int, help='The number of samples.', default = 2000 )
@click.option( '--nvariants', type=int, help='The number of variants.', default = 10000 )
@click.option( '--nancestors', type=int, help='The number of ancestors (low number means high relatedness, high number means low relatedness, default = 1000).', default = 1000 )
@click.option( '--nsegments', type=int, help='Average number of independently inherited segments (high number means low LD, low number means high LD).' )
@click.option( '--create-pair/--no-create-pair', help='Create a .pair file in the output prefix that contains all possible pairs of variants.', default = False )
@click.option( '--out', help='Output plink file.', type=click.Path( writable = True ), required = True )
def epigen(maf, nsamples, nvariants, nancestors, nsegments, create_pair, out):
    generate.write_related( nvariants, nsamples, nancestors, nsegments, out, maf = maf, create_pair = create_pair )
