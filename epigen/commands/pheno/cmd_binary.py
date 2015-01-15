import click

from epigen.commands.command import CommandWithHelp

@click.command( 'binary', cls = CommandWithHelp, short_help='Generates binary phenotypes for given plink data.' )
def epigen():
    pass
