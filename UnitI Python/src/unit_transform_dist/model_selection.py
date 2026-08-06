"""
model_selection.py
===================

Automatic selection of the best-fitting base distribution for a given
(0,1) sample: fits every candidate distribution (from the package registry,
or a user-supplied list) via MLE and ranks them by AIC/BIC/log-likelihood.
"""

from typing import List, Optional
import logging
import numpy as np
import numpy.typing as npt

from .core import UnitTransformDist, fit_unit_transform
from .extras import aic, bic, log_likelihood, ks_test
from .generate_all import get_positive_distributions

logger = logging.getLogger(__name__)

__all__ = ["auto_select_distribution"]


def auto_select_distribution(data: npt.ArrayLike, candidates: Optional[List[str]] = None,
                              criterion: str = "aic", top_n: Optional[int] = None):
    """
    Fit every candidate base distribution to `data` via MLE and rank the
    results by information criterion.

    Parameters
    ----------
    data : array_like
        Observations in (0,1).
    candidates : list of str, optional
        Names of scipy.stats positive-support distributions to try.
        Defaults to every such distribution known to the package registry.
    criterion : {'aic', 'bic', 'loglik'}, default 'aic'
        Ranking criterion (ascending for aic/bic, descending for loglik).
    top_n : int, optional
        If given, only return the top N results.

    Returns
    -------
    results : list of dict
        Each dict has keys: 'name', 'dist', 'aic', 'bic', 'loglik',
        'ks_stat', 'ks_pvalue'. Sorted best-first.
    """
    data = np.asarray(data, dtype=float)
    if np.any((data <= 0) | (data >= 1)):
        raise ValueError("All data must lie in the open interval (0,1).")

    if candidates is None:
        candidates = get_positive_distributions()

    results = []
    for name in candidates:
        try:
            fitted = fit_unit_transform(data, name)
            n_params = fitted.base_dist.numargs + 2
            ll = log_likelihood(fitted, data)
            aic_val = aic(fitted, data, n_params)
            bic_val = bic(fitted, data, n_params)
            try:
                ks_stat, ks_p = ks_test(fitted, data)
            except Exception:
                ks_stat, ks_p = np.nan, np.nan
            if not np.isfinite(ll):
                continue
            results.append({
                "name": name,
                "dist": fitted,
                "loglik": ll,
                "aic": aic_val,
                "bic": bic_val,
                "ks_stat": ks_stat,
                "ks_pvalue": ks_p,
            })
        except Exception as e:
            logger.debug(f"Skipping {name}: {e}")
            continue

    if not results:
        raise RuntimeError("No candidate distribution could be fit to the data.")

    reverse = criterion == "loglik"
    key = {"aic": "aic", "bic": "bic", "loglik": "loglik"}.get(criterion)
    if key is None:
        raise ValueError("criterion must be one of 'aic', 'bic', 'loglik'.")
    results.sort(key=lambda r: r[key], reverse=reverse)

    if top_n is not None:
        results = results[:top_n]
    return results
