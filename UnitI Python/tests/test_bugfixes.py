import numpy as np
import pytest
from unit_transform_dist import make_unit_transform, fit_moments, UnitTransformExtras, fit_unit_transform


def test_fit_moments_weibull_min_does_not_crash():
    """Regression test: fit_moments('weibull_min', ...) previously called
    scipy.stats.gamma(...) (a distribution object) instead of the gamma
    function, raising a TypeError."""
    true = make_unit_transform("weibull_min", base_shape_params=(2.0,))
    data = true.rvs(size=300, random_state=42)
    fitted = fit_moments(data, "weibull_min")
    assert fitted.base_scale > 0


def test_summary_prints_and_returns_string(capsys):
    """Regression test: summary() previously used logger.info() with no
    handler configured, so nothing was printed and None was returned."""
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=200, random_state=42)
    fitted = fit_unit_transform(data, "gamma")
    model = UnitTransformExtras(fitted, data)
    report = model.summary()
    captured = capsys.readouterr()
    assert report is not None
    assert "Summary" in captured.out
    assert "AIC" in report


def test_summary_silent_when_verbose_false(capsys):
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=100, random_state=1)
    fitted = fit_unit_transform(data, "gamma")
    model = UnitTransformExtras(fitted, data)
    report = model.summary(verbose=False)
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "AIC" in report


def test_profile_likelihood_ci_no_nameerror():
    """Regression test: profile_likelihood_ci previously called
    scipy.optimize.brentq without importing the `scipy` top-level module."""
    from unit_transform_dist import profile_likelihood_ci
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=300, random_state=42)
    fitted = fit_unit_transform(data, "gamma")
    lower, upper = profile_likelihood_ci(fitted, data, param_index=0)
    assert lower != lower or lower < upper  # allow nan, but no crash


def test_profile_likelihood_ci_correct_threshold_and_nuisance_profiling():
    """Regression test for two compounding bugs:
    (1) chdtri(v, p) args were passed as chdtri(1-alpha, 1) instead of
        chdtri(1, alpha), collapsing the chi-square threshold to ~0.
    (2) nuisance parameters (loc, scale) were fixed at their full-MLE
        values instead of being re-optimised for each candidate value of
        the parameter of interest, understating the true profile-likelihood
        interval width.
    Both together caused the 95% CI to degenerate to (MLE, MLE) and to
    frequently exclude the true parameter value.
    """
    from unit_transform_dist import profile_likelihood_ci
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=500, random_state=42)
    lower, upper = profile_likelihood_ci(true, data, param_index=0, alpha=0.05)
    assert np.isfinite(lower) and np.isfinite(upper)
    assert upper - lower > 0.1  # interval must not collapse to a point
    assert lower <= 2.5 <= upper
