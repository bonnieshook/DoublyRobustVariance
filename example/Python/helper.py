#####################################################################################################
# Helper functions
#
# Paul Zivich (2024/01/24)
#####################################################################################################

import numpy as np


def probability_to_odds(prob):
    """Helper function to convert a probability into an odds

    Parameters
    ----------
    prob : float, ndarray
    """
    return prob / (1 - prob)


def bound_unit(val, mini, maxi):
    """Helper function to bound a variable between [0,1].

    Parameters
    ----------
    val : ndarray, int, float
        Value(s) to bound.
    mini : float, int
        Minimum value for the bounding procedure
    maxi : float, int
        Maximum value for the bounding procedure
    """
    return (val - mini) / (maxi - mini)


def unbound_unit(val, mini, maxi):
    """Helper function to unbound a variable that is bound [0,1] using some range.

    Parameters
    ----------
    val : ndarray, int, float
        Value(s) to bound.
    mini : float, int
        Minimum value for the unbounding procedure
    maxi : float, int
        Maximum value for the unbounding procedure
    """
    return val * (maxi - mini) + mini
