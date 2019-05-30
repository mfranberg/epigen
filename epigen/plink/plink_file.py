from plinkio import plinkfile

class PlinkFile:
    ##
    # Constructor.
    #
    # @param path Prefix to the plink file.
    # @param phenotype Phenotypes of all individuals.
    # @param is_binary Determines whether phenotype should be interpreted
    #                  as binary or not.
    # @param iid_prefix Prefix for iids.
    #
    def __init__(self, path, phenotype, is_binary = True, iid_prefix = "iid"):
        samples = [ ]
        for i, p in enumerate( phenotype ):
            iid = "fid{0}".format( i )
            fid = "{0}{1}".format( iid_prefix, i )

            if is_binary:
                sample = plinkfile.Sample( iid, fid, "0", "0", 1, p )
                samples.append( sample )
            else:
                sample = plinkfile.Sample( iid, fid, "0", "0", 1, 3, p )
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

