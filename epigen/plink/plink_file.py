from plinkio import plinkfile

class PlinkFile:
    ##
    # Constructor.
    #
    # @param path Prefix to the plink file.
    # @param ncases Number of cases.
    # @param ncontrols Number of controls.
    #
    def __init__(self, path, ncases, ncontrols):
        samples = [ ]
        for i in range( ncases + ncontrols ):
            phenotype = int( i < ncases + 1 )

            iid = "iid{0}".format( i )
            fid = "fid{0}".format( i )
            affection = int( i < ncases + 1 )

            sample = plinkfile.Sample( iid, fid, "0", "0", 1, affection )
            samples.append( sample )

        self.plink_file = plinkfile.create( path, samples )
        self.index = 1

    ##
    # Writes the given row to the plink file.
    #
    # @param i Index of the variant.
    # @param row The genotypes.
    #
    def write(self, i, row):
        name = "rs{0}".format( i )
        bp_position = i
        locus = plinkfile.Locus( 1, name, 0, i, "A", "G" )

        self.plink_file.write_row( locus, row )

    ##
    # Closes the plink file.
    #
    def close(self):
        self.plink_file.close( )
        pass

