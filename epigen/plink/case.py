##
# Helper class that writes a list of pairs
# and whether they are an interaction or not.
#
class CaseFile:
    ##
    # Constructor.
    #
    # @param path Path to the file.
    #
    def __init__(self, path):
        self.file = open( path, "w" )

    ##
    # Writes a pair with the next indices and info on
    # whether they are an interaction or not.
    #
    # @param is_case Boolean that indicates whether this is a real interaction
    #                or not.
    #
    def write(self, pair, is_case):
        self.file.write( "{0} {1}\n".format( pair, int( is_case ) ) )

    ##
    # Closes the file.
    #
    def close(self):
        self.file.close( )
