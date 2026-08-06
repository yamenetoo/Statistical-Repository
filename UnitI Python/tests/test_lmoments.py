import numpy as np
import pytest
from unit_transform_dist import make_unit_transform, sample_lmoments, lmom_ratios, fit_lmoments


def test_sample_lmoments_basic():
    rng = np.random.default_rng(0)
    x = rng.uniform(0.01, 0.99, size=500)
    l = sample_lmoments(x, n_moments=4)
    assert len(l) == 4
    assert np.isclose(l[0], np.mean(x), atol=1e-8)
    assert l[1] > 0  # L-scale must be positive


def test_lmom_ratios_keys():
    rng = np.random.default_rng(1)
    x = rng.beta(2, 5, size=300)
    ratios = lmom_ratios(x)
    assert set(ratios.keys()) == {"l1", "l2", "t3", "t4"}


def test_fit_lmoments_recovers_reasonable_params():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=800, random_state=42)
    fitted = fit_lmoments(data, "gamma")
    assert fitted.base_shape_params[0] > 0
    # sanity: fitted pdf integrates close to 1
    y = np.linspace(0.001, 0.999, 2000)
    integral = np.trapezoid(fitted.pdf(y), y)
    assert abs(integral - 1.0) < 0.05
