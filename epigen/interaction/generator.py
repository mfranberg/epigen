import itertools

##
# Class that generates all possible interactions under
# a given null model. It does so by generating all possible
# marginal vectors (p1, p2, p3) at two loci, and  combines
# them with the given null function to yield a complete
# penetrance (p1, p2, p3, p4, p5, p6, p7, p8, p9).
#
class InteractionGenerator:
    ##
    # @param mat_add The null model to use.
    #
    def __init__(self, mat_add):
        self.mat_add = mat_add

    ##
    # Converts a marginal null vector to a
    # complete penetrance matrix for two loci.
    #
    # e.g.
    #
    # (p1, p2, p3) -> ( p1, p2, p3, p1, p2, p3, p1, p2, p3 )
    #
    # or
    #
    # (p1, p2, p3) -> ( p1, p1, p1, p2, p2, p2, p3, p3, p3 )
    #
    # @param a The marginal vector. 
    # @param other If false each row has the same penetrance,
    #              otherwise each column.
    #
    # @return The penetrance matrix as a vector.
    #
    def to_mat(self, a, other = False):
        l = list( )
        for i in range(3):
            for j in range(3):
                index = i
                if other:
                    index = j

                l.append( a[ index ] )

        return tuple( l )

    ##
    # Generates all possible interaction and null models by
    # generating all possible marginal vectors and combining
    # them with the given null function, to see which null 
    # models can be reached. The interactions are then defined
    # as all models minus these null models.
    #
    def generate(self):
        all_mats = set( a for a in itertools.product( (0,1), repeat = 9 ) )
        null_mats = set( )
        for a in itertools.product( (0,1), repeat = 3 ):
            for b in itertools.product( (0,1), repeat = 3 ):
                a_mat = self.to_mat( a, True )
                b_mat = self.to_mat( b, True )
                null_mats.add( self.mat_add( a_mat, b_mat ) )

                a_mat = self.to_mat( a, True )
                b_mat = self.to_mat( b, False )
                null_mats.add( self.mat_add( a_mat, b_mat ) )
                
                a_mat = self.to_mat( a, False )
                b_mat = self.to_mat( b, True )
                null_mats.add( self.mat_add( a_mat, b_mat ) )

                a_mat = self.to_mat( a, False )
                b_mat = self.to_mat( b, False )
                null_mats.add( self.mat_add( a_mat, b_mat ) )

        int_mats = all_mats.difference( null_mats )

        return ( int_mats, null_mats )


##
# The OR null model, if either locus is one the disease
# develops.
#
def mat_or(a, b):
    return tuple( map( lambda x: x[ 0 ] | x[ 1 ], zip (a,b) ) )
