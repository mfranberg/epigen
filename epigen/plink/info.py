from epigen.plink import variant

import json

##
# Computes the heritability V(P|G) / V(P).
#
# @param model The type of model used (normal, binomial, poisson etc)
# @param mu The mean value for each genotype.
# @param maf The minor allele frequency for both locus.
# @param dispersion The dispersion parameter if applicable.
#
# @return The heritability.
#
def compute_heritability(model, mu, maf, dispersion):
    p = [ ( 1 - maf[ 0 ] )**2, 2 * maf[ 0 ] * ( 1 - maf[ 0 ] ), ( maf[ 0 ] )**2 ]
    q = [ ( 1 - maf[ 1 ] )**2, 2 * maf[ 1 ] * ( 1 - maf[ 1 ] ), ( maf[ 1 ] )**2 ]

    joint_maf =  [ p[ 0 ] * q[ 0 ], p[ 0 ] * q[ 1 ], p[ 0 ] * q[ 2 ],
                   p[ 1 ] * q[ 0 ], p[ 1 ] * q[ 1 ], p[ 1 ] * q[ 2 ],
                   p[ 2 ] * q[ 0 ], p[ 2 ] * q[ 1 ], p[ 2 ] * q[ 2 ] ]

    pop_mu = sum( m * f for m, f in zip( mu, joint_maf ) )
    h_numerator = sum( f * m**2 for m, f in zip( mu, joint_maf ) ) - pop_mu**2

    h_denominator = h_numerator
    if model == "binomial":
        h_denominator += sum( m*(1-m) * f for m, f in zip( mu, joint_maf ) )
    elif model == "poisson":
        h_denominator += pop_mu
    else:
        h_denominator += sum( d**2 * f for d, f in zip( [ dispersion ] * 9, joint_maf ) )

    H2 = h_numerator / h_denominator
    if abs( H2 ) <= 1e-6:
        return 0.0
    else:
        return H2

##
# Computes the disease prevalence E[P]
#
# @param mu The mean value for each genotype.
# @param maf The minor allele frequency for both variants.
#
# @return The disease prevalence.
#
def compute_prevalence(mu, maf):
    p = [ ( 1 - maf[ 0 ] )**2, 2 * maf[ 0 ] * ( 1 - maf[ 0 ] ), ( maf[ 0 ] )**2 ]
    q = [ ( 1 - maf[ 1 ] )**2, 2 * maf[ 1 ] * ( 1 - maf[ 1 ] ), ( maf[ 1 ] )**2 ]

    joint_maf =  [ p[ 0 ] * q[ 0 ], p[ 0 ] * q[ 1 ], p[ 0 ] * q[ 2 ],
                   p[ 1 ] * q[ 0 ], p[ 1 ] * q[ 1 ], p[ 1 ] * q[ 2 ],
                   p[ 2 ] * q[ 0 ], p[ 2 ] * q[ 1 ], p[ 2 ] * q[ 2 ] ]

    return sum( m * f for m, f in zip( mu, joint_maf ) )

##
# Make a monte-carlo estimate of the prevalence E[P]
#
# @param mu_map The mapping from a specific genotype to a mean value.
# @param mafs The minor allele frequency for all variants.
#
# @return The estimated disease prevalence.
#
def estimate_prevalence(mu_map, mafs):
    num_iter = 1000
    mean_value = 0.0
    for i in range( num_iter ):
        variants = variant.generate_variant_set( mafs )
        mean_value += mu_map.map( variants ) / num_iter

    return mean_value

##
# Make a monte-carlo estimate of the heritability.
#
# @param model The type of model used (normal, binomial, poisson etc)
# @param mu The mean value for each genotype.
# @param maf The minor allele frequency for both locus.
# @param dispersion The dispersion parameter if applicable.
# @param pop_mu The population mean.
#
def estimate_heritability(model, mu_map, maf, dispersion, pop_mu):
    pop_var = 0.0
    if model == "binomial":
        pop_var = pop_mu * ( 1 - pop_mu )
    elif model == "normal":
        pop_var = dispersion
    elif model == "poisson":
        pop_var = pop_mu

    num_iter = 1000
    mu_list = list( )
    mean_gen_var = 0.0
    for i in range( num_iter ):
        variants = variant.generate_variant_set( maf )
        mu = mu_map.map( variants )
        mu_list.append( mu )
        mean_gen_var += mu / num_iter

    gen_var = (1.0/(num_iter - 1)) * sum( (m - mean_gen_var)**2 for m in mu_list )

    return gen_var / ( gen_var + pop_var )

##
# Writes some information about the model to a json file.
#
# @param model The type of model (binomial, normal, etc)
# @param mu The mean value for each genotype, or a mapping from
#           a specific genotype to a specific mean value.
# @param maf The minor allele frequency.
# @param dispersion The dispersion.
# @param sample_size The sample size.
# @param output_path The output path of the json file.
# @param info Additional parameters to write.
#
def write_info(model, mu, maf, dispersion, sample_size, output_path, info = dict( )):
    if model == "binomial":
        info[ "case-control-ratio" ] = float( sample_size[ 0 ] ) / sum( sample_size )

    if len( maf ) == 2:
        info[ "prevalence" ] = compute_prevalence( mu, maf )
        info[ "heritability" ] = compute_heritability( model, mu, maf, dispersion )
    else:
        pop_mu = estimate_prevalence( mu, maf )
        info[ "prevalence" ] = pop_mu
        info[ "heritability" ] = estimate_heritability( model, mu, maf, dispersion, pop_mu )

    with open( output_path, "w" ) as info_file:
        json.dump( info, info_file )
