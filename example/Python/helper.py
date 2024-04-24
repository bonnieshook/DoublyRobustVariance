#####################################################################################################
# Helper functions
#
# Paul Zivich (2024/01/24)
#####################################################################################################

import numpy as np
from scipy.stats import norm
from delicatessen.utilities import inverse_logit, logit


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


def inf_func_inference(theta, y, a, PSM, OM, OM1, OM0, unbound_y=False):
    """Function to compute the influence-function-based variance given the estimated parameters and same inputs as the
    the estimating equations.

    Parameters
    ----------
    theta : ndarray, list, vector
        *Solved* parameter vector
    y : ndarray, list, vector
        1-dimensional vector of n observed values.
    a : ndarray, list, vector
        1-dimensional vector of n observed values. The A values should all be 0 or 1.
    PSM : ndarray, list, vector
        2-dimensional vector of n observed values for the propensity score model.
    OM : ndarray, list, vector
        2-dimensional vector of n observed values for the outcome model.
    OM1 : ndarray, list, vector
        2-dimensional vector of n observed values for the outcome model under the action plan where ``a=1``.
    OM0 : ndarray, list, vector, None, optional
        2-dimensional vector of n observed values for the outcome model under the action plan where ``a=0``.
    unbound_y : bool, optional
        Whether the outcome predictions need to be unbounded. Default is False. Should only be set to True for TMLE.
    """
    # Breaking parameters into subsets
    id_ps = PSM.shape[1] + 3           # Number of parameters in PS model
    id_om = OM.shape[1] + id_ps        # Number of parameters in outcome model
    ace, mu1, mu0 = theta[:3]          # Parameters of interest
    alpha = theta[3:id_ps]             # PS model parameters
    beta = theta[id_ps:id_om]          # Outcome model parameters

    # Generating predictions
    pi = inverse_logit(np.dot(PSM, alpha))   # Generate PS
    y1hat = np.dot(OM1, beta)                # Generate Y | A=1 prediction
    y0hat = np.dot(OM0, beta)                # Generate Y | A=0 prediction

    # Optional unbounding step for TMLE
    if unbound_y:
        continuous_bound = 1e-5                            # Note: assumes default from ee_tmle
        min_y = np.min(y) - continuous_bound               # Get lower value
        max_y = np.max(y) + continuous_bound               # Get upper value
        y1p = inverse_logit(theta[-2] + logit(y1hat))      # Targeted Y1
        y0p = inverse_logit(theta[-1] + logit(y0hat))      # Targeted Y0
        y1hat = unbound_unit(y1p, mini=min_y, maxi=max_y)  # Unbound Y1
        y0hat = unbound_unit(y0p, mini=min_y, maxi=max_y)  # Unbound Y0

    # Computing the influence-function-based variance
    y1_tilde = (y*a/pi - y1hat*(a-pi)/pi)                  # Y1 contribution
    y0_tilde = (y*(1-a)/(1-pi) + y0hat*(a-pi)/(1-pi))      # Y0 contribution
    inf_func = y1_tilde - y0_tilde - ace                   # Influence function
    var = np.var(inf_func, ddof=1) / y.shape[0]            # Variance of influence function
    se = var ** 0.5                                        # Convert to SE

    # Computing the confidence intervals
    z_alpha = norm.ppf(1 - 0.05 / 2, loc=0, scale=1)       # Critical value for 95% CI
    lcl = ace - z_alpha * se                               # Lower 95% CI
    ucl = ace + z_alpha * se                               # Upper 95% CI

    # Return standard error and 95% CI
    return se, lcl, ucl
