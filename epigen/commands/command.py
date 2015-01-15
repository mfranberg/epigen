import click
import os
import sys
from epigen import commands

class ComplexCLI( click.MultiCommand ):
    def __init__(self, cmd_subdir, *args, **kwargs):
        self.cmd_subdir = cmd_subdir
        self.cmd_dir = os.path.join( os.path.abspath( os.path.dirname( commands.__file__ ) ), cmd_subdir )
        super( ComplexCLI, self ).__init__( *args, **kwargs )

    def list_commands(self, ctx):
        rv = []
        for filename in os.listdir( self.cmd_dir ):
            if filename.startswith( "cmd_" ) and filename.endswith( ".py" ):
                rv.append( self.cmd_subdir + "-" + filename[ 4:-3 ] )

        rv.sort( )

        return rv

    def get_command(self, ctx, name):
        try:
            if sys.version_info[ 0 ] == 2:
                name = name.encode( "ascii", "replace" )

            name = name.replace( self.cmd_subdir + "-", "" )
            mod_name = "epigen.commands.{0}.cmd_{1}".format( self.cmd_subdir, name )
            mod = __import__( mod_name, None, None, [ "epigen" ] )
        except ImportError:
            return

        return mod.epigen

class CommandWithHelp(click.Command):
    """A Command subclass that adds the help automatically.
    """

    def parse_args(self, ctx, args):
        """Parse arguments sent to this command.
        The code for this method is taken from MultiCommand:
        https://github.com/mitsuhiko/click/blob/master/click/core.py
        It is Copyright (c) 2014 by Armin Ronacher.
        See the license:
        https://github.com/mitsuhiko/click/blob/master/LICENSE
        """
        if not args and not ctx.resilient_parsing:
            click.echo( ctx.get_help( ) )
            ctx.exit( )
        return super( CommandWithHelp, self ).parse_args( ctx, args )
