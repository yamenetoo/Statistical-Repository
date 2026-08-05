import numpy as np
import pytest
from unit_transform_dist import (
    make_unit_transform, fit_unit_transform, UnitTransformDist,
    pearson_residuals, cox_snell_residuals, quantile_residuals, cooks_distance,
)
from scipy import stats


def test_pearson_residuals_shape_and_finiteness():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=200, random_state=42)
    fitted = fit_unit_transform(data, "gamma")
    r = pearson_residuals(fitted, data)
    assert r.shape == data.shape
    assert np.all(np.isfinite(r))


def test_cox_snell_residuals_approx_exponential():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=2000, random_state=1)
    r = cox_snell_residuals(true, data)
    # Cox-Snell residuals under the true model ~ Exponential(1): mean ~ 1
    assert 0.7 < np.mean(r) < 1.3


def test_quantile_residuals_approx_normal():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=3000, random_state=2)
    r = quantile_residuals(true, data)
    assert abs(np.mean(r)) < 0.15
    assert abs(np.std(r) - 1.0) < 0.15


def test_cooks_distance_flags_outlier():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=60, random_state=3)
    data = np.clip(data, 1e-4, 1 - 1e-4)
    data_with_outlier = data.copy()
    data_with_outlier[0] = 0.999  # extreme value
    infl = cooks_distance(UnitTransformDist, stats.gamma, data_with_outlier)
    assert infl.shape == data_with_outlier.shape
    # the injected outlier should be among the most influential points
    assert np.nanargmax(infl) == 0
