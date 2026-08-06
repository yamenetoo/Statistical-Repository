"""
generate_all.py
================
Automatically builds a transformed version (Y = 1/(1+X)) for every positive‑support
(0, inf) distribution available in scipy.stats, and collects them into a single
dictionary ready to be used as a package of statistical distributions on (0,1).
"""

import numpy as np
from scipy import stats
from .core import UnitTransformDist


def get_positive_distributions():
    """Returns the names of all continuous distributions in scipy.stats whose support is (0, inf)."""
    names = []
    for name in dir(stats):
        obj = getattr(stats, name)
        if isinstance(obj, stats.rv_continuous):
            try:
                a, b = obj.a, obj.b
            except Exception:
                continue
            if a == 0 and b == np.inf:
                names.append(name)
    return sorted(names)


def build_all_unit_transforms(default_shape_value=1.0):
    """
    Builds a dictionary {dist_name: UnitTransformDist object} for all positive‑support
    distributions. Uses a simple default value (1.0) for every required shape parameter —
    the user can later replace the parameters by manually constructing a UnitTransformDist
    with different parameters.

    Note: Some distributions (e.g., studentized_range) require special constraints on
    parameters and may fail with default values — they are automatically skipped and
    recorded in failed[].
    """
    distributions = {}
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
            # Quick sanity check to ensure the methods actually work
            _ = dist.pdf(0.5)
            _ = dist.cdf(0.5)
            _ = dist.ppf(0.5)
            distributions[name] = dist
        except Exception as e:
            failed[name] = str(e)

    return distributions, failed


# Pre‑built dictionary that can be loaded when the module is imported
# (optional – you may comment this out if building all distributions at once
# is undesirable for performance reasons)
ALL_UNIT_DISTRIBUTIONS = None
FAILED_DISTRIBUTIONS = None


def load_all(force=False):
    """Builds (or reuses) the dictionary of all transformed distributions."""
    global ALL_UNIT_DISTRIBUTIONS, FAILED_DISTRIBUTIONS
    if ALL_UNIT_DISTRIBUTIONS is None or force:
        ALL_UNIT_DISTRIBUTIONS, FAILED_DISTRIBUTIONS = build_all_unit_transforms()
    return ALL_UNIT_DISTRIBUTIONS, FAILED_DISTRIBUTIONS