##
# Helper class that writes a list of pairs
# and the index of the model that they are
# generated from.
#
class ModelFile:
    ##
    # Constructor.
    #
    # @param path Path to the file.
    #
    def __init__(self, path):
        self.file = open( path, "w" )

    ##
    # Writes a pair with the next indices and the index
    # of the model.
    #
    # @param snp1 Index of snp1.
    # @param snp2 Index of snp2.
    # @param model_index Index of the model.
    #
    def write(self, pair, model_index):
        self.file.write( "{0} {1}\n".format( pair, int( model_index ) ) )

    ##
    # Closes the file.
    #
    def close(self):
        self.file.close( )
