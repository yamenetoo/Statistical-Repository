"""
censored.py
===========

Maximum-likelihood fitting of UnitTransformDist under right-, left-, or
interval-censoring, and under truncation of the observation window.
Censoring/truncation are common in survival- and proportion-type data
(e.g. a fraction reported only as "<0.05" or "at least 0.9").
"""

from typing import Optional, Tuple
import numpy as np
import numpy.typing as npt
from scipy import stats
from scipy.optimize import minimize

from .core import UnitTransformDist, FittingError

__all__ = ["fit_censored"]

# Censoring status codes (Turnbull-style convention)
EXACT = 0
RIGHT = 1   # known only to exceed the recorded value
LEFT = 2    # known only to be below the recorded value
INTERVAL = 3  # known to lie within [value, value2]


def fit_censored(data: npt.ArrayLike, base_dist_name: str,
                  status: Optional[npt.ArrayLike] = None,
                  upper: Optional[npt.ArrayLike] = None,
                  truncation: Optional[Tuple[float, float]] = None,
                  x0: Optional[npt.ArrayLike] = None) -> UnitTransformDist:
    """
    Fit a UnitTransformDist to censored and/or truncated data by MLE.

    Parameters
    ----------
    data : array_like
        Recorded values in (0,1). For interval-censored observations this is
        the lower bound of the interval.
    base_dist_name : str
        Name of a scipy.stats positive-support distribution.
    status : array_like of int, optional
        Censoring code per observation: 0=exact, 1=right-censored,
        2=left-censored, 3=interval-censored. Defaults to all-exact.
    upper : array_like, optional
        Upper bound of the interval, required where status==3 (interval).
    truncation : (float, float), optional
        If given, the observation window (a, b) subset of (0,1): the
        likelihood is normalised by [F(b) - F(a)] to account for the fact
        that only observations within (a, b) could have been recorded
        (truncation, as opposed to censoring).
    x0 : array_like, optional
        Initial parameter guess (shape_params..., loc, scale).

    Returns
    -------
    UnitTransformDist
        The fitted distribution.
    """
    data = np.asarray(data, dtype=float)
    n = len(data)
    if status is None:
        status = np.zeros(n, dtype=int)
    else:
        status = np.asarray(status, dtype=int)
    if upper is None:
        upper = np.full(n, np.nan)
    else:
        upper = np.asarray(upper, dtype=float)

    base_dist = getattr(stats, base_dist_name)
    n_shape = base_dist.numargs
    n_params = n_shape + 2

    finite_mask = status != RIGHT  # exact/left/interval lower bounds are informative
    x_est = 1.0 / np.clip(data[finite_mask], 1e-6, 1 - 1e-6) - 1.0
    if x0 is None:
        x0 = np.ones(n_params)
        x0[-2] = 0.0
        x0[-1] = max(np.std(x_est), 1e-3) if len(x_est) else 1.0

    def neg_log_likelihood(params: np.ndarray) -> float:
        shape_params = tuple(params[:n_shape])
        loc, scale = params[-2], params[-1]
        if scale <= 0:
            return 1e10
        try:
            dist = UnitTransformDist(base_dist=base_dist, base_shape_params=shape_params,
                                      base_loc=loc, base_scale=scale, name="temp_cens")
            ll = 0.0

            exact_mask = status == EXACT
            if np.any(exact_mask):
                lp = dist.logpdf(data[exact_mask])
                if not np.all(np.isfinite(lp)):
                    return 1e10
                ll += np.sum(lp)

            right_mask = status == RIGHT
            if np.any(right_mask):
                sf = dist.sf(data[right_mask])
                sf = np.clip(sf, 1e-300, 1.0)
                ll += np.sum(np.log(sf))

            left_mask = status == LEFT
            if np.any(left_mask):
                cdf = dist.cdf(data[left_mask])
                cdf = np.clip(cdf, 1e-300, 1.0)
                ll += np.sum(np.log(cdf))

            interval_mask = status == INTERVAL
            if np.any(interval_mask):
                lo = data[interval_mask]
                hi = upper[interval_mask]
                if np.any(np.isnan(hi)):
                    return 1e10
                prob = dist.cdf(hi) - dist.cdf(lo)
                prob = np.clip(prob, 1e-300, 1.0)
                ll += np.sum(np.log(prob))

            if truncation is not None:
                a, b = truncation
                trunc_prob = dist.cdf(b) - dist.cdf(a)
                trunc_prob = max(trunc_prob, 1e-300)
                ll -= n * np.log(trunc_prob)

            return -ll if np.isfinite(ll) else 1e10
        except Exception:
            return 1e10

    loc_lower = -0.999 * np.min(x_est) if len(x_est) else -1.0
    bounds = [(1e-8, None)] * n_shape + [(loc_lower, None), (1e-8, None)]
    result = minimize(neg_log_likelihood, x0, bounds=bounds, method='L-BFGS-B')
    if not result.success:
        # Nelder-Mead fallback: more robust to poor scaling / flat regions,
        # at the cost of speed and a formal convergence guarantee.
        result = minimize(neg_log_likelihood, x0, method='Nelder-Mead',
                           options={"maxiter": 5000, "xatol": 1e-8, "fatol": 1e-8})
        if not result.success:
            raise FittingError(f"Censored/truncated MLE failed: {result.message}")

    params = result.x
    shape_params = tuple(params[:n_shape])
    loc, scale = params[-2], params[-1]
    return UnitTransformDist(base_dist=base_dist, base_shape_params=shape_params,
                              base_loc=loc, base_scale=scale,
                              name=f"unit_{base_dist_name}_censored")
