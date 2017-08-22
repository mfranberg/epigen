import click
from plinkio import plinkfile
import random

from epigen.plink import generate, genmodels, info
from epigen.plink.util import find_rows, sample_loci_set, find_beta0, generate_beta, compute_mafs
from epigen.commands.command import CommandWithHelp

@click.command( 'multiple', cls = CommandWithHelp, short_help='Generates environmental variables.' )
@click.argument( 'plink_file', type=click.Path( ) )
@click.option( '--num-variables', type=int, help='The number of environmental variables to generate.', default = 1 )
@click.option( '--out', type = click.File( 'w' ), help='Output phenotype file.', required = True )
def epigen(plink_file, num_variables, out):
    input_file = plinkfile.open( plink_file ) 
    samples = [ (s.fid, s.iid) for s in  input_file.get_samples( ) ]

    generate.write_environment( samples, num_variables, out )
