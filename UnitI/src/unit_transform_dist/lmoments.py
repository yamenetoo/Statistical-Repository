"""
lmoments.py
===========

L-moments computation (sample L-moments) and L-moment based fitting for
UnitTransformDist objects. L-moments are linear combinations of order
statistics; they are more robust than conventional moments and are the
natural choice for bounded-support distributions such as those produced
by the unit transform Y = 1/(1+X).

References
----------
Hosking, J.R.M. (1990). "L-moments: Analysis and Estimation of Distributions
using Linear Combinations of Order Statistics." JRSS B, 52(1), 105-124.
"""

from typing import Tuple, Optional
import numpy as np
import numpy.typing as npt
from scipy.optimize import minimize

from .core import UnitTransformDist, FittingError

__all__ = ["sample_lmoments", "lmom_ratios", "fit_lmoments"]


def sample_lmoments(data: npt.ArrayLike, n_moments: int = 4) -> np.ndarray:
    """
    Compute the first `n_moments` sample L-moments (l1, l2, l3, l4, ...)
    using the unbiased probability-weighted-moments estimator.

    Parameters
    ----------
    data : array_like
        1-D sample.
    n_moments : int, default 4
        Number of L-moments to compute (>=1).

    Returns
    -------
    lmoments : ndarray, shape (n_moments,)
    """
    x = np.sort(np.asarray(data, dtype=float))
    n = len(x)
    if n < n_moments:
        raise ValueError("Sample size must be >= n_moments.")

    # Unbiased probability weighted moments b_r = (1/n) * sum_i C(i-1,r)/C(n-1,r) * x_(i)
    b = np.zeros(n_moments)
    i = np.arange(1, n + 1)
    b[0] = np.mean(x)
    for r in range(1, n_moments):
        num = np.ones(n)
        den = 1.0
        for k in range(r):
            num = num * (i - 1 - k)
            den *= (n - 1 - k)
        num = np.where(num < 0, 0.0, num)
        b[r] = np.sum(num * x) / (n * den)

    # Convert PWMs to L-moments using shifted Legendre coefficients
    l = np.zeros(n_moments)
    l[0] = b[0]
    if n_moments > 1:
        l[1] = 2 * b[1] - b[0]
    if n_moments > 2:
        l[2] = 6 * b[2] - 6 * b[1] + b[0]
    if n_moments > 3:
        l[3] = 20 * b[3] - 30 * b[2] + 12 * b[1] - b[0]
    return l


def lmom_ratios(data: npt.ArrayLike, n_moments: int = 4) -> dict:
    """
    Returns a dict with l1 (location/mean), l2 (scale), t3 (L-skewness),
    t4 (L-kurtosis).
    """
    l = sample_lmoments(data, n_moments=max(n_moments, 4))
    out = {"l1": l[0], "l2": l[1]}
    out["t3"] = l[2] / l[1] if l[1] != 0 else np.nan
    out["t4"] = l[3] / l[1] if l[1] != 0 else np.nan
    return out


def _theoretical_lmoments(dist: UnitTransformDist, n_moments: int = 4,
                           n_quad: int = 256) -> np.ndarray:
    """
    Approximate theoretical L-moments of `dist` via Gauss-Legendre quadrature
    on the quantile function: L_r = integral_0^1 Q(u) P*_{r-1}(u) du,
    where P*_{r-1} are shifted Legendre polynomials.
    """
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    # map nodes from [-1,1] to [0,1]
    u = 0.5 * (nodes + 1.0)
    w = 0.5 * weights
    q = dist.ppf(u)

    def shifted_legendre(r, u):
        # P*_r(u) = P_r(2u - 1)
        return np.polynomial.legendre.Legendre.basis(r)(2 * u - 1)

    l = np.zeros(n_moments)
    for r in range(n_moments):
        l[r] = np.sum(w * q * shifted_legendre(r, u))
    return l


def fit_lmoments(data: npt.ArrayLike, base_dist_name: str,
                  n_moments: int = 4, x0: Optional[npt.ArrayLike] = None,
                  n_quad: int = 128) -> UnitTransformDist:
    """
    Fit a UnitTransformDist by matching sample L-moments to theoretical
    L-moments (computed numerically via the quantile function), minimising
    the squared distance between the two vectors. Works for any base
    distribution, unlike Method of Moments which only supports a handful.

    Parameters
    ----------
    data : array_like
        Observations in (0, 1).
    base_dist_name : str
        Name of a scipy.stats continuous distribution with positive support.
    n_moments : int, default 4
        Number of L-moments to match (2, 3, or 4 recommended).
    x0 : array_like, optional
        Initial guess for (shape_params..., loc, scale). Defaults to a
        reasonable heuristic.
    n_quad : int
        Number of quadrature nodes used for theoretical L-moments.

    Returns
    -------
    UnitTransformDist
        Fitted distribution.
    """
    from scipy import stats

    data = np.asarray(data, dtype=float)
    if np.any((data <= 0) | (data >= 1)):
        raise ValueError("All data must be in the open interval (0,1).")

    base_dist = getattr(stats, base_dist_name)
    n_shape = base_dist.numargs
    n_params = n_shape + 2
    n_moments = max(n_moments, n_params)

    target = sample_lmoments(data, n_moments=n_moments)

    if x0 is None:
        x0 = np.ones(n_params)
        x0[-2] = 0.0   # loc
        x0[-1] = 1.0   # scale

    def objective(params):
        shape_params = tuple(params[:n_shape])
        loc, scale = params[-2], params[-1]
        if scale <= 0:
            return 1e10
        try:
            trial = UnitTransformDist(base_dist=base_dist, base_shape_params=shape_params,
                                       base_loc=loc, base_scale=scale, name="temp_lmom")
            theo = _theoretical_lmoments(trial, n_moments=n_moments, n_quad=n_quad)
            return float(np.sum((theo - target) ** 2))
        except Exception:
            return 1e10

    bounds = [(1e-6, None)] * n_shape + [(None, None), (1e-6, None)]
    result = minimize(objective, x0, bounds=bounds, method='L-BFGS-B')
    if not result.success:
        raise FittingError(f"L-moment optimisation failed: {result.message}")

    params = result.x
    shape_params = tuple(params[:n_shape])
    loc, scale = params[-2], params[-1]
    return UnitTransformDist(base_dist=base_dist, base_shape_params=shape_params,
                              base_loc=loc, base_scale=scale,
                              name=f"unit_{base_dist_name}_lmom")
