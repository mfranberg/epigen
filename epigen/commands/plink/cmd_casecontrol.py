import click
from plinkio import plinkfile
import random
import json

from epigen.util import probability
from epigen.plink import generate, genmodels, info
from epigen.plink.util import generate_beta
from epigen.commands.command import CommandWithHelp

def generate_mafs(maf, n):
    generate_maf = lambda: random.betavariate( 0.4679562, 0.4679562 )
    if maf:
        generate_maf = lambda: maf[ 0 ] + ( maf[ 1 ] - maf[ 0 ] ) * random.random( )

    return [ generate_maf( ) for i in range( n ) ]

@click.command( 'additive', cls = CommandWithHelp, short_help='Generate case/control data with both true and false variants.' )
@click.option( '--maf', nargs=2, type=probability.probability, help='If set MAF is generated uniformly between these two values (default use exp distribution).', default = None )
@click.option( '--mu', nargs=9, type=float, help='Space-separated list of floating point numbers that represents the mean value for each genotype, specified row-wise from left to right.', default = None )
@click.option( '--beta0', type=float, help='Sets the intercept.', default = 0.0 )
@click.option( '--beta-sim', nargs=2, type=float, help='The mean and variance of the beta variables (taken from a normal).', default = None )
@click.option( '--beta', nargs=9, type=float, help='Space-separated list of regression coefficients a, b1, b2, g1, g2, d11, d12, d21 and d22.', default = None )
@click.option( '--link', type=click.Choice( genmodels.get_links( ).keys( ) ), help="The link function to use.", default = "default" )
@click.option( '--dispersion', type=float, help="The dispersion parameter to use.", default=1.0 )
@click.option( '--num-true', type=int, help='The number of loci that is involved in the phenotype (used in --beta-sim).', default = 2 )
@click.option( '--num-false', type=int, help='The number of loci that is not involved in the phenotype', default = 10 )
@click.option( '--sample-size', nargs=2, type=int, help='Number of samples (if only one group only first argument will be used).', default = [2000, 2000] )
@click.option( '--out', type = click.Path( exists = False ), help='Output prefix (pheno will be .pheno).', required = True )
def epigen(maf, mu, beta0, beta, beta_sim, link, dispersion, num_true, num_false, sample_size, out): 
    pheno_generator = None

    if (mu or beta) and num_true != 2:
        print( "epigen: error: With 'mu' or 'beta' --num-true must be 2." )
        exit( 1 )

    mu_values = mu
    if mu and not beta and not beta_sim:
        mu_map = genmodels.GeneralMuMap( mu )
        pheno_generator = genmodels.get_pheno_generator( "binomial", mu_map, dispersion )
    elif not mu and beta and not beta_sim:
        lf = genmodels.get_link( "binomial", link )
        mu_values = genmodels.get_mean_values( beta, lf )
    
        mu_map = genmodels.GeneralMuMap( mu_values )
        pheno_generator = genmodels.get_pheno_generator( "binomial", mu_map, dispersion )
    elif not mu and not beta and beta_sim:
        gen_beta = generate_beta( num_true, beta_sim[ 0 ], beta_sim[ 1 ] )
        mu_values = genmodels.AdditiveMuMap( beta0, gen_beta, genmodels.get_link( "binomial", link ) )
        pheno_generator = genmodels.get_pheno_generator( "binomial", mu_values, dispersion )
    else:
        print( "epigen: error: Only one of --mu, --beta-sim or --beta must be set." )
        exit( 1 )
    
    mafs = generate_mafs( maf, num_true + num_false )
    info.write_info( "binomial", mu_values, mafs[ :num_true ], dispersion, sample_size, out + ".info", { "num-true" : num_true, "num-false" : num_false } )
    with open( out + ".pheno", "w" ) as pheno_file:
        generate.write_casecontrol_data( pheno_generator, sample_size, mafs, num_true, num_false, out, pheno_file, False ) 


