##
# Helper class that writes a list of snp pairs.
#
class PairFile:
    ##
    # Constructor.
    #
    # @param path Path to the file.
    #
    def __init__(self, path):
        self.file = open( path, "w" )

    ##
    # Writes a pair with the next index.
    #
    def write(self, pair):
        self.file.write( pair + "\n" )

    ##
    # Closes the file.
    #
    def close(self):
        self.file.close( )
