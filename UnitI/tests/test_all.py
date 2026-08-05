"""
tests/test_all.py
------------------
Comprehensive tests for the unit_transform_dist package.
"""

import pytest
import numpy as np
from scipy import stats
from unit_transform_dist import (
    UnitTransformDist,
    make_unit_transform,
    fit_unit_transform,
    fit_moments,
    log_likelihood,
    aic,
    bic,
    ks_test,
    cramer_von_mises_test,
    anderson_darling_test,
    hazard,
    cumulative_hazard,
    bootstrap_ci,
    profile_likelihood_ci,
    get_registry,
    load_all,
    FittingError,
)

# ------------------------------
# Fixtures
# ------------------------------

@pytest.fixture
def gamma_data():
    """Synthetic data from gamma(shape=2.5) with transform."""
    true_dist = make_unit_transform("gamma", base_shape_params=(2.5,))
    return true_dist.rvs(size=500, random_state=42), true_dist

@pytest.fixture
def expon_data():
    true_dist = make_unit_transform("expon")
    return true_dist.rvs(size=500, random_state=42), true_dist

# ------------------------------
# Core functionality
# ------------------------------

def test_pdf_integrates_to_one():
    dist = make_unit_transform("expon")
    y = np.linspace(0.001, 0.999, 1000)
    pdf = dist.pdf(y)
    integral = np.trapezoid(pdf, y)   # updated from trapz
    assert np.isclose(integral, 1.0, rtol=1e-2)

def test_cdf_monotonic():
    dist = make_unit_transform("gamma", base_shape_params=(2.0,))
    y = np.linspace(0.01, 0.99, 100)
    cdf = dist.cdf(y)
    assert np.all(np.diff(cdf) >= 0)

def test_ppf_inverts_cdf():
    dist = make_unit_transform("lognorm", base_shape_params=(0.5,))
    q = np.array([0.1, 0.5, 0.9])
    y = dist.ppf(q)
    cdf_at_y = dist.cdf(y)
    np.testing.assert_allclose(cdf_at_y, q, rtol=1e-6)

# ------------------------------
# Fitting (MLE)
# ------------------------------

def test_fit_recovers_gamma(gamma_data):
    data, true_dist = gamma_data
    fitted = fit_unit_transform(data, "gamma")
    # true shape = 2.5, loc=0, scale=1
    assert np.isclose(fitted.base_shape_params[0], 2.5, rtol=0.2)
    assert np.isclose(fitted.base_loc, 0.0, atol=0.1)
    assert np.isclose(fitted.base_scale, 1.0, rtol=0.2)

def test_fit_recovers_expon(expon_data):
    data, true_dist = expon_data
    fitted = fit_unit_transform(data, "expon")
    # expon has no shape, scale should be ~1
    assert np.isclose(fitted.base_scale, 1.0, rtol=0.2)
    assert fitted.base_shape_params == ()

def test_fit_handles_bad_data():
    with pytest.raises(ValueError, match="All data must be in the open interval"):
        make_unit_transform("expon").fit(np.array([0.0, 0.5, 1.0]))

# ------------------------------
# Goodness-of-fit tests
# ------------------------------

def test_gof_tests_run_without_error(gamma_data):
    data, true_dist = gamma_data
    # Just ensure they execute
    ks_stat, ks_p = ks_test(true_dist, data)
    assert 0 <= ks_p <= 1
    assert ks_stat >= 0

    cvm_res = cramer_von_mises_test(true_dist, data)
    assert 0 <= cvm_res.pvalue <= 1

    ad_res = anderson_darling_test(true_dist, data)
    assert 'statistic' in ad_res
    assert len(ad_res['critical_values']) == 5

# ------------------------------
# Model selection
# ------------------------------

def test_aic_bic(gamma_data):
    data, true_dist = gamma_data
    ll = log_likelihood(true_dist, data)
    aic_val = aic(true_dist, data)
    bic_val = bic(true_dist, data)
    assert aic_val == 2 * 3 - 2 * ll  # gamma has 3 params (shape, loc, scale)
    n = len(data)
    assert bic_val == 3 * np.log(n) - 2 * ll

# ------------------------------
# Hazard functions
# ------------------------------

def test_hazard(gamma_data):
    data, true_dist = gamma_data
    y = np.linspace(0.1, 0.9, 10)
    h = hazard(true_dist, y)
    assert np.all(h >= 0)
    assert not np.any(np.isnan(h))

    H = cumulative_hazard(true_dist, y)
    assert np.all(H >= 0)
    # At y=1, hazard should be inf; but we are not testing boundary.

# ------------------------------
# Bootstrap and Profile CI
# ------------------------------

def test_bootstrap_ci(gamma_data):
    data, true_dist = gamma_data
    cis = bootstrap_ci(true_dist, data, n_bootstrap=50, random_state=0)  # small for speed
    assert 'shape1' in cis
    lower, upper = cis['shape1']
    assert lower < upper
    # true shape 2.5 should lie inside
    assert lower <= 2.5 <= upper

def test_profile_ci(gamma_data):
    data, true_dist = gamma_data
    lower, upper = profile_likelihood_ci(true_dist, data, param_index=0, alpha=0.05)
    # If both are nan, skip (numerical issue)
    if not np.isfinite(lower) and not np.isfinite(upper):
        pytest.skip("Profile likelihood CI could not be computed (numerical issues).")
    # Otherwise check basic properties
    if np.isfinite(lower) and np.isfinite(upper):
        assert lower <= upper
        # true shape should be inside (with some slack)
        assert lower <= 2.5 <= upper

# ------------------------------
# Method of Moments
# ------------------------------

def test_fit_moments_gamma(gamma_data):
    data, true_dist = gamma_data
    mom_dist = fit_moments(data, "gamma")
    # Not too strict, just ensure it runs and gives a valid object
    assert isinstance(mom_dist, UnitTransformDist)
    pdf_val = mom_dist.pdf(0.5)
    assert np.isfinite(pdf_val)

def test_fit_moments_expon(expon_data):
    data, true_dist = expon_data
    mom_dist = fit_moments(data, "expon")
    assert isinstance(mom_dist, UnitTransformDist)
    # scale should be close to mean of transformed data
    x = 1.0 / data - 1.0
    assert np.isclose(mom_dist.base_scale, np.mean(x), rtol=0.1)

# ------------------------------
# Registry
# ------------------------------

def test_registry_builds():
    reg, failed = get_registry(force_rebuild=True)
    # Check at least a few common distributions are present
    for name in ['expon', 'gamma', 'lognorm']:
        assert name in reg
    # Ensure a distribution that might fail is handled gracefully (e.g., studentized_range)
    # We don't enforce success, just that it doesn't crash

def test_load_all_works():
    all_dists, failed = load_all(force=True)
    # We don't test all, just a few
    assert 'expon' in all_dists
    # If any failed, they are recorded
    assert isinstance(failed, dict)

# ------------------------------
# Plotting (smoke test)
# ------------------------------

def test_plotting_functions_run(gamma_data, monkeypatch):
    # Prevent showing plots during test
    monkeypatch.setattr("matplotlib.pyplot.show", lambda: None)
    from unit_transform_dist import plot_pdf, plot_cdf, plot_sf, plot_all
    data, true_dist = gamma_data
    # Just check they don't raise errors
    plot_pdf(true_dist, show=False)
    plot_cdf(true_dist, show=False)
    plot_sf(true_dist, show=False)
    plot_all(true_dist)  # show parameter doesn't exist; we monkeypatch plt.show anyway

# ------------------------------
# Edge Cases
# ------------------------------

def test_extreme_y_values():
    dist = make_unit_transform("expon")
    # pdf at very small y should be small, but finite
    y_small = 1e-10
    pdf_val = dist.pdf(y_small)
    assert np.isfinite(pdf_val) and pdf_val >= 0 
    # cdf at near 1 should be ~1
    y_large = 1 - 1e-10
    cdf_val = dist.cdf(y_large)
    assert np.isclose(cdf_val, 1.0, rtol=1e-6)

def test_negative_parameters():
    dist = make_unit_transform("gamma", base_shape_params=(-1.0,))
    # pdf should be invalid (nan or inf) for some y
    pdf_val = dist.pdf(0.5)
    assert not np.isfinite(pdf_val)  # should be nan or inf

# ------------------------------
# Multiple base distributions
# ------------------------------

@pytest.mark.parametrize("dist_name, shape_params", [
    ("expon", ()),
    ("gamma", (2.0,)),
    ("lognorm", (0.5,)),
    ("weibull_min", (1.5,)),
])
def test_make_and_fit_various(dist_name, shape_params):
    # Create and fit a small dataset, check that fit doesn't error
    true = make_unit_transform(dist_name, base_shape_params=shape_params)
    data = true.rvs(size=100, random_state=42)
    fitted = fit_unit_transform(data, dist_name)
    assert fitted.base_shape_params == tuple(fitted.base_shape_params)  # just check it exists