import click
from epigen.commands.command import CommandWithHelp
from epigen.plink import generate

@click.command( 'data', cls = CommandWithHelp, short_help="Generates a plink file without a phenotype." )
@click.option( '--maf', nargs=2, type=generate.probability, help='If set MAF is generated uniformly between these two values (default use exp distribution).', default = None )
@click.option( '--nsamples', type=int, help='The number of samples.', default = 2000 )
@click.option( '--nvariants', type=int, help='The number of variants.', default = 10000 )
@click.option( '--out', help='Output plink file.', type=click.Path( writable = True ), required = True )
def epigen(maf, nsamples, nvariants, out):
    generate.write_single( nvariants, nsamples, out, maf = maf )
