"""
diagnostics.py
==============

Residual analysis and influence diagnostics for fitted UnitTransformDist
models: Pearson and (Cox-Snell) deviance-type residuals, randomised
quantile residuals, and an approximate case-deletion (Cook's-distance-like)
influence measure based on the empirical influence function.
"""

from typing import Optional
import numpy as np
import numpy.typing as npt
from scipy import stats as spstats

from .core import UnitTransformDist

__all__ = [
    "pearson_residuals", "cox_snell_residuals", "quantile_residuals",
    "cooks_distance",
]


def pearson_residuals(dist: UnitTransformDist, data: npt.ArrayLike) -> np.ndarray:
    """(y_i - E[Y]) / sd(Y), using the fitted distribution's overall mean/std.

    Note: this is the classical Pearson residual for a homogeneous model
    (constant parameters). For regression models with observation-specific
    scale, compute residuals directly from `UnitTransformRegression`.
    """
    data = np.asarray(data, dtype=float)
    mu = dist.mean()
    sigma = dist.std()
    if sigma == 0 or not np.isfinite(sigma):
        raise ValueError("Distribution has zero or undefined variance.")
    return (data - mu) / sigma


def cox_snell_residuals(dist: UnitTransformDist, data: npt.ArrayLike) -> np.ndarray:
    """
    Cox-Snell residuals: r_i = -log(S(y_i)) = cumulative hazard at y_i.
    Under a correctly specified model these should behave like a standard
    Exponential(1) sample.
    """
    data = np.asarray(data, dtype=float)
    sf = np.clip(dist.sf(data), 1e-300, 1.0)
    return -np.log(sf)


def quantile_residuals(dist: UnitTransformDist, data: npt.ArrayLike,
                        random_state: Optional[int] = None) -> np.ndarray:
    """
    Randomised (Dunn-Smyth) quantile residuals: r_i = Phi^{-1}(F(y_i)),
    which are approximately standard Normal under a correctly specified
    continuous model. Useful for QQ-plots against N(0,1).
    """
    data = np.asarray(data, dtype=float)
    u = np.clip(dist.cdf(data), 1e-12, 1 - 1e-12)
    return spstats.norm.ppf(u)


def cooks_distance(dist_class, base_dist, data: npt.ArrayLike,
                    fit_kwargs: Optional[dict] = None) -> np.ndarray:
    """
    Approximate Cook's-distance-like influence measure via case deletion:
    for each observation i, refit the MLE without it and measure the
    (Euclidean, standardised) shift in fitted parameters relative to the
    full-data fit.

    This is O(n) refits and therefore intended for moderate sample sizes
    (roughly n <= a few thousand).

    Parameters
    ----------
    dist_class : class
        Typically `UnitTransformDist`.
    base_dist : scipy.stats.rv_continuous instance
        e.g. `scipy.stats.gamma`.
    data : array_like
        Observations in (0,1).
    fit_kwargs : dict, optional
        Passed through to `.fit()`.

    Returns
    -------
    ndarray, shape (n,)
        Influence score per observation (larger = more influential).
    """
    data = np.asarray(data, dtype=float)
    n = len(data)
    fit_kwargs = fit_kwargs or {}

    full = dist_class(base_dist=base_dist, name="full")
    full_params = np.array(full.fit(data, **fit_kwargs))

    influence = np.zeros(n)
    for i in range(n):
        subset = np.delete(data, i)
        try:
            temp = dist_class(base_dist=base_dist, name="loo")
            params_i = np.array(temp.fit(subset, **fit_kwargs))
            diff = params_i - full_params
            scale = np.where(np.abs(full_params) > 1e-8, np.abs(full_params), 1.0)
            influence[i] = np.sqrt(np.sum((diff / scale) ** 2))
        except Exception:
            influence[i] = np.nan
    return influence
