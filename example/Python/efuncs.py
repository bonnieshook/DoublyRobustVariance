#####################################################################################################
# Estimating functions for each of the three doubly robust estimators
#
# Paul Zivich (2024/01/24)
#####################################################################################################

import numpy as np
from delicatessen.estimating_equations import ee_regression, ee_aipw
from delicatessen.utilities import logit, inverse_logit

from helper import bound_unit, unbound_unit


def ee_aipw_plugin(theta, y, a, PSM, OM, OM1, OM0):
    """Plug-in AIPW estimator based on the built-in functionality from delicatessen.

    Parameters
    ----------
    theta : ndarray, list, vector
        Parameter vector
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
    """
    # Plug-in AIPW: uses the built-in estimatng equation from delicatessen
    return ee_aipw(theta=theta, y=y, A=a, W=PSM, X=OM, X1=OM1, X0=OM0)


def ee_aipw_wreg(theta, y, a, PSM, OM, OM1, OM0):
    """Weighted regression AIPW estimator. Estimating function is coded from scratch, but uses built-in regression
    functionalities from delicatessen. Note that this code only works for continuous outcomes.

    Parameters
    ----------
    theta : ndarray, list, vector
        Parameter vector
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
    """
    # Divide parameters into their respective models
    id_ps = PSM.shape[1] + 3         # Number of parameters in nuisance PS model
    ace, mu1, mu0 = theta[0:3]       # Parameters of interest
    alpha = theta[3:id_ps]           # PS model parameters
    beta = theta[id_ps:]             # Outcome model parameters

    # Propensity score model
    ee_act = ee_regression(theta=alpha,         # Logistic model for the PS
                           X=PSM, y=a,
                           model='logistic')
    pi = inverse_logit(np.dot(PSM, alpha))      # Computing the PS
    ipw = a / pi + (1-a) / (1-pi)               # Compute the IPW from PS

    # IPW-Outcome model
    ee_out = ee_regression(theta=beta,          # OLS for the outcome
                           X=OM, y=y,
                           model='linear',
                           weights=ipw)         # ... with IPW as weights
    y1hat = np.dot(OM1, beta)                   # Predicted Y under a=1
    y0hat = np.dot(OM0, beta)                   # Predicted Y under a=0

    # Weighted regression AIPW estimator
    ee_mu1 = y1hat - mu1                               # EE mu1
    ee_mu0 = y0hat - mu0                               # EE mu0
    ee_ace = np.ones(y.shape[0]) * (mu1 - mu0) - ace   # EE average causal effect

    # Return stacked estimating equations
    return np.vstack([ee_ace, ee_mu1, ee_mu0, ee_act, ee_out])


def ee_tmle(theta, y, a, PSM, OM, OM1, OM0, continuous_bound=1e-5):
    """TMLE. Estimating function is coded from scratch, but uses built-in regression functionalities from delicatessen.
    Note that this code only works for continuous outcomes.

    Parameters
    ----------
    theta : ndarray, list, vector
        Parameter vector
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
    continuous_bound : float, optional
        Values to shift the lower and upper bounds by.
    """
    # Divide parameters into their respective models
    id_ps = PSM.shape[1] + 3        # Number of parameters in nuisance PS model
    id_om = OM.shape[1] + id_ps     # Number of parameters in nuisance outcome model
    ace, mu1, mu0 = theta[0:3]      # Parameters of interest
    alpha = theta[3:id_ps]          # PS model parameters
    beta = theta[id_ps:id_om]       # Outcome model parameters
    eta1 = [theta[id_om], ]         # Y1 targeting model parameter
    eta0 = [theta[id_om+1], ]       # Y0 targeting model parameter

    # Some meta parameters for the bounding and unbounding
    min_y = np.min(y) - continuous_bound  # Shifted minimum
    max_y = np.max(y) + continuous_bound  # Shifted maximum

    # Bounding outcomes and saving bounding info
    y = bound_unit(y, mini=min_y, maxi=max_y)

    # Propensity score model
    ee_act = ee_regression(theta=alpha,           # Logistic model for PS
                           X=PSM, y=a,
                           model='logistic')
    pi = inverse_logit(np.dot(PSM, alpha))        # Computing the PS
    clever1 = a / pi                              # Clever covariate for A=1
    clever0 = (1-a) / (1-pi)                      # Clever covariate for A=0

    # Outcome model
    ee_out = ee_regression(theta=beta,            # OLS for outcome
                           X=OM, y=y,
                           model='linear')
    y1hat = np.dot(OM1, beta)                     # Predicted Y under a=1
    y0hat = np.dot(OM0, beta)                     # Predicted Y under a=0

    # Targeting models
    intercept = np.ones([y.shape[0], 1])          # Intercept-only design matrix
    ee_trg1 = ee_regression(theta=eta1,           # Targeting model for Y1
                            X=intercept, y=y,     # ... intercept only for Y
                            offset=logit(y1hat),  # ... offset of prediction
                            weights=clever1,      # ... clever cov weight
                            model='logistic')
    ee_trg0 = ee_regression(theta=eta0,           # Targeting model for Y0
                            X=intercept, y=y,     # ... intercept only for Y
                            offset=logit(y0hat),  # ... offset of prediction
                            weights=clever0,      # ... clever cov weight
                            model='logistic')

    # Applying TMLE for estimates
    y1tilde = inverse_logit(eta1 + logit(y1hat))     # Targeted Y1 prediction
    y0tilde = inverse_logit(eta0 + logit(y0hat))     # Targeted Y0 prediction

    # Parameters of interest
    ee_mu1 = unbound_unit(y1tilde, mini=min_y, maxi=max_y) - mu1  # EE mu1
    ee_mu0 = unbound_unit(y0tilde, mini=min_y, maxi=max_y) - mu0  # EE mu0
    ee_ace = np.ones(y.shape)*(mu1 - mu0) - ace                   # EE average causal effect

    # Return stacked estimating equations
    return np.vstack([ee_ace, ee_mu1, ee_mu0, ee_act, ee_out, ee_trg1, ee_trg0])
