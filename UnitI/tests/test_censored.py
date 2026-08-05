import numpy as np
import pytest
from unit_transform_dist import make_unit_transform, fit_censored


def test_fit_censored_all_exact_matches_uncensored():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=400, random_state=42)
    fitted = fit_censored(data, "gamma")  # default: all exact
    assert 1.5 < fitted.base_shape_params[0] < 4.0


def test_fit_censored_with_right_censoring():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=400, random_state=1)
    status = np.zeros(len(data), dtype=int)
    # right-censor the smallest 10% of values (interpreted as "known to exceed")
    threshold_idx = np.argsort(data)[:40]
    status[threshold_idx] = 1
    fitted = fit_censored(data, "gamma", status=status)
    assert fitted.base_shape_params[0] > 0


def test_fit_censored_with_truncation():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=2000, random_state=2)
    mask = (data > 0.1) & (data < 0.9)
    trunc_data = data[mask]
    fitted = fit_censored(trunc_data, "gamma", truncation=(0.1, 0.9))
    assert fitted.base_shape_params[0] > 0
