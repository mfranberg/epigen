##
# Converts the given value to a floating point value
# if it represents a probability, otherwise it raises
# a ValueError.
#
# @param prob_string String representing a probability.
#
# @return String as a floating point probability.
#
def probability(prob_string):
    prob_value = float( prob_string )
    if 0.0 <= prob_value <= 1.0:
        return prob_value
    else:
        raise ValueError( "Probability not between 0 and 1." )

    return prob_value
