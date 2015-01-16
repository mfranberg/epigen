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
# Writes the normal phenotypes for two snps and a set of individuals to a file.
#
# @param sample_list List of plinkio.plinkfile.Sample.
# @param snp1_row List of genotypes for the first snp.
# @param snp2_row List of genotypes for the second snp.
# @param model Mean matrix as a length 9 list
# @param sd Common standard deviation
# @param output_file The phenotypes will be written to this file.
#
def write_phenotypes(sample_list, snp1_row, snp2_row, model, sd, output_file, plink_format = False):
    output_file.write( "FID\tIID\tPheno\n" )

    num_samples = 0
    for sample, snp1, snp2 in zip( sample_list, snp1_row, snp2_row ):
        pheno = generate.generate_phenotype_cont( snp1, snp2, model, sd )
        pheno_str = "NA"
        if plink_format:
            pheno_str = "-9"

        if pheno != None:
            pheno_str = str( pheno )
            num_samples += 1

        output_file.write( "{0}\t{1}\t{2}\n".format( sample.fid, sample.iid, pheno_str ) )

    return num_samples

@click.command( 'cont', cls = CommandWithHelp, short_help='Generates binary phenotypes for given plink data.' )
@click.option( '--model', nargs=9, type=generate.probability, help='Space-separated list of floating point numbers that represents the mean for each genotype, specified row-wise from left to right.' )
@click.option( '--sd', metavar='sd', type=float, help='Common standard deviation', required = True )
@click.option( '--pair', nargs=2, type=str, help='Name of two SNPs for which the phenotype should be based on (otherwise random).', default = None )
@click.option( "--plink-format/--no-plink-format", help="Use plink format for the phenotype file.", default = False )
@click.option( '--out', type=click.Path(writable = True), help='Output phenotype file.', required = True )
@click.argument( 'plink_file', type=click.Path( ) )
def epigen(plink_file, model, sd, pair, plink_format, out):
    input_file = plinkfile.open( plink_file )
    snp1_index, snp2_index = generate.sample_two_loci( input_file.get_loci( ) )
    if pair:
        snp1_index, snp2_index = find_pair( input_file.get_loci( ), args.pair[ 0 ], args.pair[ 1 ] )

    snp1_row, snp2_row = find_genotypes( input_file, snp1_index, snp2_index )

    with open( out, "w" ) as output_file:
        num_samples = write_phenotypes( input_file.get_samples( ), snp1_row, snp2_row, model, sd, output_file, plink_format = plink_format )
        print "Wrote {0} samples".format( num_samples ) 
