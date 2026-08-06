"""
bayes.py
========

Lightweight Bayesian estimation for UnitTransformDist: MAP estimation and a
random-walk Metropolis-Hastings sampler. No external MCMC dependency is
required (only numpy/scipy), so this stays usable in restricted environments.
"""

from typing import Callable, Dict, Optional, Tuple
import numpy as np
import numpy.typing as npt
from scipy.optimize import minimize

from .core import UnitTransformDist, FittingError

__all__ = ["map_estimate", "metropolis_hastings", "PriorSpec"]


class PriorSpec:
    """
    Simple container mapping parameter index -> log-density function.

    Example
    -------
    >>> priors = PriorSpec()
    >>> priors.add(0, lambda v: -0.5 * ((v - 1.0) / 2.0) ** 2)  # weak Normal(1, 2) prior on shape
    """

    def __init__(self):
        self._log_priors: Dict[int, Callable[[float], float]] = {}

    def add(self, index: int, log_density: Callable[[float], float]) -> "PriorSpec":
        self._log_priors[index] = log_density
        return self

    def log_prior(self, params: npt.ArrayLike) -> float:
        total = 0.0
        for idx, fn in self._log_priors.items():
            total += fn(params[idx])
        return total


def _log_posterior(params, base_dist, data, n_shape, priors: Optional[PriorSpec]):
    shape_params = tuple(params[:n_shape])
    loc, scale = params[-2], params[-1]
    if scale <= 0 or any(p <= 0 for p in shape_params) if n_shape else scale <= 0:
        return -np.inf
    try:
        dist = UnitTransformDist(base_dist=base_dist, base_shape_params=shape_params,
                                  base_loc=loc, base_scale=scale, name="temp_bayes")
        ll = np.sum(dist.logpdf(data))
        if not np.isfinite(ll):
            return -np.inf
        lp = priors.log_prior(params) if priors is not None else 0.0
        return ll + lp
    except Exception:
        return -np.inf


def map_estimate(data: npt.ArrayLike, base_dist_name: str,
                  priors: Optional[PriorSpec] = None,
                  x0: Optional[npt.ArrayLike] = None) -> UnitTransformDist:
    """
    Maximum a posteriori estimate: like MLE, but maximises log-likelihood
    plus log-prior density. With `priors=None` this reduces to MLE.
    """
    from scipy import stats
    data = np.asarray(data, dtype=float)
    base_dist = getattr(stats, base_dist_name)
    n_shape = base_dist.numargs
    n_params = n_shape + 2

    if x0 is None:
        x0 = np.ones(n_params)
        x0[-2] = 0.0
        x0[-1] = 1.0

    def neg_log_post(params):
        lp = _log_posterior(params, base_dist, data, n_shape, priors)
        return -lp if np.isfinite(lp) else 1e10

    bounds = [(1e-8, None)] * n_shape + [(None, None), (1e-8, None)]
    result = minimize(neg_log_post, x0, bounds=bounds, method='L-BFGS-B')
    if not result.success:
        raise FittingError(f"MAP optimisation failed: {result.message}")

    params = result.x
    shape_params = tuple(params[:n_shape])
    loc, scale = params[-2], params[-1]
    return UnitTransformDist(base_dist=base_dist, base_shape_params=shape_params,
                              base_loc=loc, base_scale=scale,
                              name=f"unit_{base_dist_name}_map")


def metropolis_hastings(data: npt.ArrayLike, base_dist_name: str,
                         n_samples: int = 5000, burn_in: int = 1000,
                         proposal_scale: float = 0.05,
                         priors: Optional[PriorSpec] = None,
                         x0: Optional[npt.ArrayLike] = None,
                         random_state: Optional[int] = None) -> Dict[str, np.ndarray]:
    """
    Random-walk Metropolis-Hastings sampler over (shape_params..., loc, scale),
    working on a log-transformed scale for positivity-constrained parameters
    (all shape params and scale; loc is sampled on its natural scale).

    Returns
    -------
    dict with keys:
        'samples' : ndarray (n_samples, n_params) post burn-in
        'acceptance_rate' : float
        'param_names' : list of str
    """
    from scipy import stats
    data = np.asarray(data, dtype=float)
    base_dist = getattr(stats, base_dist_name)
    n_shape = base_dist.numargs
    n_params = n_shape + 2
    rng = np.random.default_rng(random_state)

    if x0 is None:
        x0 = np.ones(n_params)
        x0[-2] = 0.0
        x0[-1] = 1.0
    current = np.array(x0, dtype=float)
    current_lp = _log_posterior(current, base_dist, data, n_shape, priors)

    positive_mask = np.array([True] * n_shape + [False, True])  # shapes & scale are positive

    total = n_samples + burn_in
    chain = np.zeros((total, n_params))
    accepted = 0

    for t in range(total):
        proposal = current.copy()
        step = rng.normal(0, proposal_scale, size=n_params)
        # random-walk on log-scale for positive parameters
        for i in range(n_params):
            if positive_mask[i]:
                proposal[i] = current[i] * np.exp(step[i])
            else:
                proposal[i] = current[i] + step[i]

        prop_lp = _log_posterior(proposal, base_dist, data, n_shape, priors)
        log_alpha = prop_lp - current_lp
        if np.log(rng.uniform()) < log_alpha:
            current, current_lp = proposal, prop_lp
            accepted += 1
        chain[t] = current

    samples = chain[burn_in:]
    param_names = [f"shape{i+1}" for i in range(n_shape)] + ["loc", "scale"]
    return {
        "samples": samples,
        "acceptance_rate": accepted / total,
        "param_names": param_names,
    }


def posterior_summary(mcmc_result: Dict[str, np.ndarray], alpha: float = 0.05) -> Dict[str, Tuple[float, float, float]]:
    """Return {param: (mean, lower, upper)} credible intervals from MCMC samples."""
    samples = mcmc_result["samples"]
    names = mcmc_result["param_names"]
    out = {}
    for i, name in enumerate(names):
        col = samples[:, i]
        out[name] = (float(np.mean(col)),
                     float(np.percentile(col, 100 * alpha / 2)),
                     float(np.percentile(col, 100 * (1 - alpha / 2))))
    return out
