import click
from plinkio import plinkfile

from epigen.plink import generate
from epigen.commands.command import CommandWithHelp

##
# Finds the indices of two specified variant names.
#
# @param loci List of locus names
# @param snp1 Snp1 name.
# @param snp2 Snp2 name.
#
# @return A tuple containing the indices of the variants.
#
def find_pair(loci, snp1, snp2):
    snp_to_index = dict( zip( map( lambda x: x.name, loci ), range( len( loci ) ) ) )
    return ( snp_to_index[ snp1 ], snp_to_index[ snp2 ] )

##
# Finds the genotypes of the two indices.
#
# @param plink_file PlinkFile object.
# @param snp1_index Index of the first variant.
# @param snp2_index Index of the second variant.
#
# @return A tuple containing the genotype vectors for
#         both variants.
#
def find_genotypes(plink_file, snp1_index, snp2_index):
    snp1_row = snp2_row = None
    for row_num, row in enumerate( plink_file ):
        if row_num == snp1_index:
            snp1_row = list( row )
        elif row_num == snp2_index:
            snp2_row = list( row )

    return ( snp1_row, snp2_row )

##
# Format the phenotype depending on whether it should
# be used with plink or not.
#
# @param pheno A binary phenotype 0, 1.
# @param plink_format True if output should be in plink format,
#                     False otherwise. Plink is (1,2,-9), non-plink
#                     is (0,1,NA).
#
# @return the phenotype in the desired format.
#
def format_phenotype(pheno, plink_format):
    if pheno in (0,1):
        if plink_format:
            return int( pheno + 1 )
        else:
            return int( pheno )
    else:
        if plink_format:
            return "-9"
        else:
            return "NA"

##
# Writes the phenotypes for two snps and a set of individuals to a file.
#
# @param sample_list List of plinkio.plinkfile.Sample.
# @param snp1_row List of genotypes for the first snp.
# @param snp2_row List of genotypes for the second snp.
# @param model Penetrance matrix as a length 9 list
# @param output_file The phenotypes will be written to this file.
# @param plink_format Use plink format for the phenotype.
#
def write_phenotypes(sample_list, snp1_row, snp2_row, model, output_file, plink_format = False):
    output_file.write( "FID\tIID\tPheno\n" )

    number_of_cases = number_of_controls = 0
    for sample, snp1, snp2 in zip( sample_list, snp1_row, snp2_row ):
        pheno = generate.generate_phenotype( snp1, snp2, model )
        if pheno == 0:
            number_of_controls += 1
        elif pheno == 1:
            number_of_cases += 1
        
        output_file.write( "{0}\t{1}\t{2}\n".format( sample.fid, sample.iid, format_phenotype( pheno, plink_format ) ) )

    return number_of_cases, number_of_controls


@click.command( 'binary', cls = CommandWithHelp, short_help='Generates binary phenotypes for given plink data.' )
@click.option( '--model', nargs=9, type=generate.probability, help='Space-separated list of floating point numbers that represents the penetrance matrix, specified row-wise from left to right.' )
@click.option( '--pair', nargs=2, type=str, help='Name of two SNPs for which the phenotype should be based on (otherwise random).', default = None )
@click.option( "--plink-format/--no-plink-format", help="Use plink format for the phenotype file.", default = False )
@click.option( '--out', type=click.Path(writable = True), help='Output phenotype file.', required = True )
@click.argument( 'plink_file', type=click.Path( ) )
def epigen(plink_file, model, pair, plink_format, out):
    input_file = plinkfile.open( plink_file ) 
    snp1_index, snp2_index = generate.sample_two_loci( input_file.get_loci( ) )
    if pair:
        snp1_index, snp2_index = find_pair( input_file.get_loci( ), args.pair[ 0 ], args.pair[ 1 ] )

    snp1_row, snp2_row = find_genotypes( input_file, snp1_index, snp2_index )

    with open( out, "w" ) as output_file:
        number_of_cases, number_of_controls = write_phenotypes( input_file.get_samples( ), snp1_row, snp2_row, model, output_file, plink_format = plink_format )
        print "Wrote {0} cases and {1} controls".format( number_of_cases, number_of_controls ) 
