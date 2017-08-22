import click
from plinkio import plinkfile
from math import sqrt
import random

from epigen.plink import generate, genmodels, info, envfile
from epigen.plink.util import find_rows, sample_loci_set, find_beta0, generate_beta, compute_mafs, sample_gxe, find_gxe, mean, stdev
from epigen.commands.command import CommandWithHelp

@click.command( 'env', cls = CommandWithHelp, short_help='Generates phenotypes using a gene-environment interaction model' )
@click.argument( 'plink_file', type=click.Path( ) )
@click.argument( 'env_file', type=click.Path( ) )
@click.option( '--beta0', type=float, help='Sets the intercept, by default it is chosen to be the mean value.', default = None )
@click.option( '--main-dist', type=float, nargs=2, help='Mean and variance for genetic main effects.', default = [0,0] )
@click.option( '--env-dist', type=float, nargs=2, help='Mean and variance for environment main effects.', default = [0,0] )
@click.option( '--gxe-dist', type=float, nargs=2, help='Mean and variance for gene-environment interaction effects.', default = [0,0] )
@click.option( '--lock-main', type=bool, help='Main effects are only generated for the interactions.', default = False )
@click.option( '--num-main', type=int, help='The number of genetic main effects (if --lock-main is set this option has no effect).', default = 1 )
@click.option( '--num-env', type=int, help='The number of environmental effects.', default = 1 )
@click.option( '--num-gxe', type=int, help='The number of gene-environment interactions.', default = 1 )
@click.option( '--model', type=click.Choice( genmodels.get_models( ) ), help="The model to use.", required = True )
@click.option( '--link', type=click.Choice( genmodels.get_links( ).keys( ) ), help="The link function to use.", default = "default" )
@click.option( '--dispersion', type=float, help="The dispersion parameter to use (if none will be remaining heritability, otherwise heritability will be rescaled).", default=None )
@click.option( '--out', type = click.File( 'w' ), help='Output phenotype file.', required=True )
def epigen(plink_file, env_file, beta0, main_dist, env_dist, gxe_dist, lock_main, num_main, num_env, num_gxe, model, link, dispersion, out):
    genotype_file = plinkfile.open( plink_file ) 
    iid = [ s.iid for s in genotype_file.get_samples( ) ]
    loci = genotype_file.get_loci( )
    snp_indices = sample_loci_set( loci, num_main )
    main_std = 0
    if num_main > 0:
        main_std = sqrt( main_dist[ 1 ] / num_main )
    genotype_beta = generate_beta( num_main, main_dist[ 0 ], main_std )

    env = envfile.openenv( env_file, iid )
    env_names = env.get_names( )
    env_indices = sample_loci_set( env_names, num_env )
    env_std = 0
    if num_env > 0:
        env_std = sqrt( env_dist[ 1 ] / num_env )
    env_beta = generate_beta( num_env, env_dist[ 0 ], env_std )

    genotype_file = plinkfile.open( plink_file ) 
    gxe_indices = sample_gxe( loci, env_names, num_gxe )
    gxe_std = 0
    if num_gxe > 0:
        gxe_std = sqrt( gxe_dist[ 1 ] / num_gxe )
    gxe_beta = generate_beta( num_gxe, gxe_dist[ 0 ], gxe_std )

    if not dispersion:
        dispersion = 1.0 - main_dist[ 1 ] - env_dist[ 1 ] - gxe_dist[ 1 ]

    gxe_data = find_gxe( genotype_file, env, snp_indices, env_indices, gxe_indices )

    all_beta = list( )
    all_beta.extend( genotype_beta )
    all_beta.extend( env_beta )
    all_beta.extend( gxe_beta )

    truth = list( )
    truth.extend( loci[ i ].name for i in snp_indices )
    truth.extend( env_names[ i ] for i in env_indices )
    truth.extend( env_names[ j ] + ":" + loci[ i ].name for i, j in gxe_indices )

    if not beta0:
        beta0 = find_beta0( gxe_data, all_beta )

    data_means = list( map( mean, gxe_data ) )
    data_stdev = list( map( stdev, gxe_data ) )

    mu_map = genmodels.AdditiveMuMap( beta0, all_beta, genmodels.get_link( model, link ), data_means, data_stdev )
    pheno_generator = genmodels.get_pheno_generator( model, mu_map, sqrt( dispersion ) )
    generate.write_general_phenotype( genotype_file.get_samples( ), gxe_data, pheno_generator, out, False )
    extra_info = { "truth" : truth }
    info.write_info( model, mu_map, None, dispersion, pheno_generator.sample_size, out.name + ".info", info = extra_info )
