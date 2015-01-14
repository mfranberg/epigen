from .plink_file import PlinkFile
from .pair import PairFile
from .case import CaseFile
from .model import ModelFile

class OutputFiles:
    def __init__(self, path, ncases, ncontrols):
        self.plink_file = PlinkFile( path, ncases, ncontrols )
        self.pair_file = PairFile( path + ".pair" )
        self.case_file = CaseFile( path + ".case" )
        self.model_file = ModelFile( path + ".model" )
        self.index = 1

    ##
    # Writes a row of genotypes to the output file.
    # 
    def write(self, row1, row2, is_case, model_index):
        self.plink_file.write( self.index, row1 )
        self.plink_file.write( self.index + 1, row2 )
        
        pair = "rs{0} rs{1}".format( self.index, self.index + 1 )
        self.pair_file.write( pair )
        self.case_file.write( pair, is_case )
        self.model_file.write( pair, model_index )

        self.index += 2

    def close(self):
        self.plink_file.close( )
        self.pair_file.close( )
        self.case_file.close( )
        self.model_file.close( )
