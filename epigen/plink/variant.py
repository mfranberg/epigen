from epigen.plink import util

##
# Generate a single variant.
#
# @param m The minor allele frequency.
#
# @return The genotype of the variant.
# 
def generate_variant(m):
    return util.sample_categorical( [ (1-m)**2, 2*m*(1-m), m**2 ], [ 0, 1, 2 ] )

##
# Generate a set of variants according to a list
# of minor allele frequencies.
#
# @param mafs A list of minor allele frequencies.
#
# @return A list of genotypes corresponding to each maf.
# 
def generate_variant_set(mafs):
    return [ generate_variant( m ) for m in mafs ]

##
# Generate a variant for a set of individuals.
#
# @param m The minor allele frequency.
# @param num_samples The number of samples.
#
# @return A genotype for each sample.
# 
def generate_variant_row(m, num_samples):
    return [ generate_variant( m ) for i in range( num_samples ) ]
