import os
import sys
import click

from epigen.commands.command import ComplexCLI

@click.command(no_args_is_help = True, cmd_subdirs = ["pair", "pheno"], cls = ComplexCLI)
def epigen():
    """Generate plink or phenotype data using an epistatic model."""

if __name__ == "__main__":
    epigen( )
