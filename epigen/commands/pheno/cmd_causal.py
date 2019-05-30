import click
from plinkio import plinkfile
import random
from math import sqrt

from epigen.plink import generate, genmodels, info
from epigen.plink.util import find_rows, sample_loci_set, find_beta0, generate_beta, compute_mafs
from epigen.commands.command import CommandWithHelp

@click.command( 'causal', cls = CommandWithHelp, short_help='Generates binary phenotypes for given plink data.' )
@click.argument( 'plink_file', type=click.Path( ) )
@click.option( '--beta0', type=float, help='Sets the intercept, by default it is chosen to get 50/50 cases and controls.', default = None )
@click.option( '--effect-h2', type=float, help='Narrow-sense heritability', required = True )
@click.option( '--effect-mean', type=float, help='Shift from zero of effect size distribution', required = True )
@click.option( '--num-causal', type=int, help='The number of loci that is involved in the phenotype.', default = 10 )
@click.option( '--model', type=click.Choice( genmodels.get_models( ) ), help="The model to use.", required = True )
@click.option( '--link', type=click.Choice( genmodels.get_links( ).keys( ) ), help="The link function to use.", default = "default" )
@click.option( '--dispersion', type=float, help="The dispersion parameter to use.", default=1.0 )
@click.option( '--out', type = click.File( 'w' ), help='Output phenotype file.', required = True )
def epigen(plink_file, beta0, effect_h2, effect_mean, num_causal, model, link, dispersion, out):
    input_file = plinkfile.open( plink_file ) 
    loci = input_file.get_loci( )
    snp_indices = sample_loci_set( loci, num_causal )
    gen_beta = generate_beta( num_causal, effect_mean, sqrt( effect_h2 / num_causal ) )
    rows = find_rows( input_file, snp_indices )

    causal_names = list( loci[ i ].name for i in snp_indices )
    name_to_beta = dict( zip( causal_names, gen_beta ) )

    # A bit weird to have 1.0 as default and change it here, but need
    # to keep it as is to avoid published scripts starting to fail
    if model == "normal" and dispersion == 1.0 and effect_h2 < 1.0:
        dispersion = sqrt( 1.0 - effect_h2 )

    if not beta0:
        beta0 = find_beta0( rows, gen_beta )

    mu_map = genmodels.AdditiveMuMap( beta0, gen_beta, genmodels.get_link( model, link ) )
    pheno_generator = genmodels.get_pheno_generator( model, mu_map, dispersion )
    generate.write_general_phenotype( input_file.get_samples( ), rows, pheno_generator, out, False )
    extra_info = { "truth" : list( loci[ i ].name for i in snp_indices ), "beta" : name_to_beta }
    info.write_info( model, mu_map, compute_mafs( rows ), dispersion, pheno_generator.sample_size, out.name + ".info", info = extra_info, multiple = True )
