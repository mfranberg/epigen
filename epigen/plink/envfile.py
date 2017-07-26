
class EnvFile:
    def __init__(self, env_path, env_file):
        self.path = env_path
        self.env_file = env_file
        self.parsed = False
        self.header = []
        self.data = []
        self.iid = dict( )

    def get_names(self):
        return self.header

    def parse(self, order):
        if self.parsed:
            return

        position = dict( (iid, i) for i, iid in  enumerate( order ) )

        header = next( self.env_file ).strip( ).split( )
        self.header = header[ 2: ]
        for i in range( len( self.header ) ):
            self.data.append( [ "NA" ] * len( order ) )

        for line in self.env_file:
            column = line.strip( ).split( )
            row = list( map( float, column[ 2: ] ) )
            iid = column[ 1 ]
            cur_pos = position[ iid ]

            for i, value in enumerate( row ):
                self.data[ i ][ cur_pos ] = value

        self.parsed = True

    def get_variables(self, indices):
        return [ self.data[ i ] for i in indices ]

    def close(self):
        self.env_file.close( )

def openenv(env_path, order):
    env_file = open( env_path, "r" )

    env = EnvFile( env_path, env_file )
    env.parse( order )

    return env
