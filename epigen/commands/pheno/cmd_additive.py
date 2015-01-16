import click
from plinkio import plinkfile
import random

from epigen.plink import generate
from epigen.plink.util import find_rows, sample_loci_set
from epigen.commands.command import CommandWithHelp

##
# Computes the second allele frequency for the
# given genotypes.
#
# @param row Genotypes 0, 1, 2. Here 3 is missing.
#
# @return The second allele frequency.
#
def compute_maf(row):
    no_missing = filter( lambda x: x != 3, row )
    p = sum( no_missing ) / ( 2.0 * len( no_missing ) )

    return p

##
# Determines a beta0 that gives the probability of being
# a case 0.5.
#
# @param rows List of genotypes for each locus.
# @param beta The beta used for each locus.
#
# @return A beta0 that makes the probability of being a case 0.5.
#
def find_beta0(rows, beta):
    maf = [ compute_maf( r ) for r in rows ]
    beta0 = -sum( b * 2 * m for b, m in zip( beta, maf ) )

    return beta0

##
# Generates effect sizes according to the effect size
# distribution.
#
# @param n The number of beta
# @param mean The mean of the beta distribution
# @param sd The standard deviation of the beta distibution
#
# @return list of effect sizes
#
def generate_beta(n, mean, sd):
    return [ random.normalvariate( mean, sd ) for i in range( n ) ]

##
# Writes the phenotypes for two snps and a set of individuals to a file.
#
# @param sample_list List of plinkio.plinkfile.Sample.
# @param snp1_row List of genotypes for the first snp.
# @param snp2_row List of genotypes for the second snp.
# @param model Penetrance matrix as a length 9 list
# @param output_file The phenotypes will be written to this file.
#
def write_phenotypes(sample_list, rows, beta, output_file, beta0 = None):
    output_file.write( "FID\tIID\tPheno\n" )

    if not beta0:
        beta0 = find_beta0( rows, beta )

    number_of_cases = number_of_controls = 0
    for i, sample in enumerate( sample_list ):
        variants = [ rows[ j ][ i ] for j in range( len( rows ) ) ]
        
        pheno = generate.generate_phenotype_additive( variants, beta0, beta )
        pheno_str = str( pheno )
        if pheno == 0:
            number_of_controls += 1
        elif pheno == 1:
            number_of_cases += 1
        else:
            pheno_str = "NA"

        output_file.write( "{0}\t{1}\t{2}\n".format( sample.fid, sample.iid, pheno_str ) )

    return number_of_cases, number_of_controls

@click.command( 'additive', cls = CommandWithHelp, short_help='Generates binary phenotypes for given plink data.' )
@click.argument( 'plink_file', type=click.Path( ) )
@click.option( '--beta0', type=float, help='Sets the intercept, by default it is chosen to get 50/50 cases and controls.', default = None )
@click.option( '--beta', nargs=2, help='The mean and variance of the beta variables (taken from a normal).', required = True )
@click.option( '--num-loci', type=int, help='The number of loci that is involved in the phenotype.', default = 10 )
@click.option( '--out', type = click.Path(writable=True), help='Output phenotype file.', required = True )
def epigen(plink_file, beta0, beta, num_loci, out):
    input_file = plinkfile.open( plink_file ) 
    loci = input_file.get_loci( )
    snp_indices = sample_loci_set( loci, num_loci )
    gen_beta = generate_beta( num_loci, beta[ 0 ], beta[ 1 ] )
    rows = find_rows( input_file, snp_indices )

    with open( out, "w" ) as output_file:
        number_of_cases, number_of_controls = write_phenotypes( input_file.get_samples( ), rows, gen_beta, output_file, beta0 = beta0 )
        print "Wrote {0} cases and {1} controls".format( number_of_cases, number_of_controls ) 
