"""
registry.py
===========

Central registry for all UnitTransformDist objects.
Provides a dictionary of distribution names → transformed distribution instances.
"""

import warnings
from scipy import stats
from .core import UnitTransformDist
from .generate_all import get_positive_distributions


# Global cache
_DISTRIBUTION_CACHE = None
_FAILED_CACHE = None
_CUSTOM_DISTS = {}


def get_registry(force_rebuild=False, default_shape_value=1.0):
    """
    Returns a dictionary {dist_name: UnitTransformDist} for all positive‑support
    distributions from scipy.stats.

    Parameters
    ----------
    force_rebuild : bool, default False
        If True, rebuild the cache from scratch.
    default_shape_value : float, default 1.0
        Default shape parameter value for distributions that require shape parameters.

    Returns
    -------
    registry : dict
        Mapping from distribution name to UnitTransformDist.
    failed : dict
        Mapping from distribution name to error message for those that failed.
    """
    global _DISTRIBUTION_CACHE, _FAILED_CACHE

    if _DISTRIBUTION_CACHE is None or force_rebuild:
        registry = {}
        failed = {}

        for name in get_positive_distributions():
            base_dist = getattr(stats, name)
            n_shapes = base_dist.numargs
            shape_params = tuple([default_shape_value] * n_shapes)

            try:
                dist = UnitTransformDist(
                    base_dist=base_dist,
                    base_shape_params=shape_params,
                    name=f"unit_{name}",
                )
                # Quick sanity check
                _ = dist.pdf(0.5)
                _ = dist.cdf(0.5)
                _ = dist.ppf(0.5)
                registry[name] = dist
            except Exception as e:
                failed[name] = str(e)

        _DISTRIBUTION_CACHE = registry
        _FAILED_CACHE = failed

    return _DISTRIBUTION_CACHE, _FAILED_CACHE


def register_custom(name, dist):
    """
    Add a custom UnitTransformDist to the registry.

    Parameters
    ----------
    name : str
        Unique name for the distribution.
    dist : UnitTransformDist
        The distribution object to store.
    """
    _CUSTOM_DISTS[name] = dist


def get_all_distributions(include_custom=True, include_failed=False):
    """
    Returns a combined dictionary of all distributions (built‑in + custom).

    Parameters
    ----------
    include_custom : bool, default True
        Include custom distributions registered via `register_custom`.
    include_failed : bool, default False
        If True, also include the failed distributions (as None) in the dictionary.

    Returns
    -------
    combined : dict
        Mapping from name to UnitTransformDist (or None if failed and included).
    """
    registry, failed = get_registry()
    combined = dict(registry)
    if include_custom:
        combined.update(_CUSTOM_DISTS)
    if include_failed:
        combined.update({name: None for name in failed})
    return combined