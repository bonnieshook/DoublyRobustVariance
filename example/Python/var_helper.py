#####################################################################################################
# Helper functions for the variance estimation
#
# Paul Zivich (2024/11/01)
#####################################################################################################

import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import norm
from formulaic import model_matrix
from delicatessen import MEstimator
from delicatessen.utilities import inverse_logit, logit

from efuncs import ee_aipw_plugin, ee_aipw_wreg, ee_tmle
from helper import bound_unit, unbound_unit, probability_to_odds


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


def aipw_plugin_bootstrap(data, y, a, a_model, y_model, iterations=200):
    """Bootstrapping procedure for plug-in AIPW"""
    data = data.copy()

    def aipw_statsmodels(d):
        # Propensity score model
        fm_a = smf.glm(a + " ~ " + a_model, d, family=sm.families.Binomial()).fit()
        pi_a = fm_a.predict()

        # Outcome model
        fm_y = smf.glm(y + " ~ " + y_model, d, family=sm.families.Gaussian()).fit()
        da = d.copy()
        da['anemia'] = 1
        y1hat = fm_y.predict(da)
        da['anemia'] = 0
        y0hat = fm_y.predict(da)

        r1 = np.mean(d[y]*d[a]/pi_a - y1hat*(d[a]-pi_a)/pi_a)
        r0 = np.mean(d[y]*(1-d[a])/(1-pi_a) + y0hat*(d[a]-pi_a)/(1-pi_a))
        return r1 - r0

    # Recomputing estimate to check everything works as intended
    ace = aipw_statsmodels(d=data)

    # Bootstrap standard error
    ests = []
    for b in range(iterations):
        ds = data.sample(n=data.shape[0], replace=True)
        estimate = aipw_statsmodels(d=ds)
        ests.append(estimate)

    se = np.std(ests, ddof=1)

    # Computing the confidence intervals
    z_alpha = norm.ppf(1 - 0.05 / 2, loc=0, scale=1)       # Critical value for 95% CI
    lcl = ace - z_alpha * se                               # Lower 95% CI
    ucl = ace + z_alpha * se                               # Upper 95% CI

    # Return standard error and 95% CI
    return se, lcl, ucl


def aipw_wreg_bootstrap(data, y, a, a_model, y_model, iterations=200):
    """Bootstrapping procedure for weighted-regression AIPW"""
    data = data.copy()

    def aipw_statsmodels(d):
        # Propensity score model
        fm_a = smf.glm(a + " ~ " + a_model, d, family=sm.families.Binomial()).fit()
        pi_a = fm_a.predict()
        iptw = d[a]/pi_a + (1-d[a])/(1-pi_a)

        # Outcome model
        fm_y = smf.glm(y + " ~ " + y_model, d, family=sm.families.Gaussian(),
                       freq_weights=iptw).fit()
        da = d.copy()
        da['anemia'] = 1
        y1hat = fm_y.predict(da)
        da['anemia'] = 0
        y0hat = fm_y.predict(da)

        r1 = np.mean(y1hat)
        r0 = np.mean(y0hat)
        return r1 - r0

    # Recomputing estimate to check everything works as intended
    ace = aipw_statsmodels(d=data)

    # Bootstrap standard error
    ests = []
    for b in range(iterations):
        ds = data.sample(n=data.shape[0], replace=True)
        estimate = aipw_statsmodels(d=ds)
        ests.append(estimate)

    se = np.std(ests, ddof=1)

    # Computing the confidence intervals
    z_alpha = norm.ppf(1 - 0.05 / 2, loc=0, scale=1)       # Critical value for 95% CI
    lcl = ace - z_alpha * se                               # Lower 95% CI
    ucl = ace + z_alpha * se                               # Upper 95% CI

    # Return standard error and 95% CI
    return se, lcl, ucl


def tmle_bootstrap(data, y, a, a_model, y_model, iterations=200):
    """Bootstrapping procedure for TMLE"""
    data = data.copy()
    continuous_bound = 1e-5
    min_y = np.min(data[y]) - continuous_bound
    max_y = np.max(data[y]) + continuous_bound
    data[y] = bound_unit(data[y], mini=min_y, maxi=max_y)

    def tmle_statsmodels(d):
        # Propensity score model
        fm_a = smf.glm(a + " ~ " + a_model, d, family=sm.families.Binomial()).fit()
        pi_a = fm_a.predict()
        clever1 = d[a] / pi_a              # Clever covariate for A=1
        clever0 = (1 - d[a]) / (1 - pi_a)  # Clever covariate for A=0

        # Outcome model
        fm_y = smf.glm(y + " ~ " + y_model, d, family=sm.families.Gaussian()).fit()
        da = d.copy()
        da['anemia'] = 1
        y1hat = fm_y.predict(da)
        da['anemia'] = 0
        y0hat = fm_y.predict(da)

        # Targeting model
        log = sm.GLM(d[y], np.repeat(1, d[y].shape[0]), offset=np.log(probability_to_odds(y1hat)), freq_weights=clever1,
                     family=sm.families.Binomial()).fit()
        eta1 = log.params[0]
        log = sm.GLM(d[y], np.repeat(1, d[y].shape[0]), offset=np.log(probability_to_odds(y0hat)), freq_weights=clever0,
                     family=sm.families.Binomial()).fit()
        eta0 = log.params[0]

        # Predicted pseudo potential outcomes
        y1tilde = inverse_logit(eta1 + logit(y1hat))  # Targeted Y1 prediction
        y1tilde = unbound_unit(y1tilde, mini=min_y, maxi=max_y)
        r1 = np.mean(y1tilde)
        y0tilde = inverse_logit(eta0 + logit(y0hat))  # Targeted Y0 prediction
        y0tilde = unbound_unit(y0tilde, mini=min_y, maxi=max_y)
        r0 = np.mean(y0tilde)
        return r1 - r0

    # Recomputing estimate to check everything works as intended
    ace = tmle_statsmodels(d=data)

    # Bootstrap standard error
    ests = []
    for b in range(iterations):
        ds = data.sample(n=data.shape[0], replace=True)
        estimate = tmle_statsmodels(d=ds)
        ests.append(estimate)

    se = np.std(ests, ddof=1)

    # Computing the confidence intervals
    z_alpha = norm.ppf(1 - 0.05 / 2, loc=0, scale=1)       # Critical value for 95% CI
    lcl = ace - z_alpha * se                               # Lower 95% CI
    ucl = ace + z_alpha * se                               # Upper 95% CI

    # Return standard error and 95% CI
    return se, lcl, ucl
