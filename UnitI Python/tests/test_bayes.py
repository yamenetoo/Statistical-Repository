import numpy as np
import pytest
from unit_transform_dist import make_unit_transform, map_estimate, metropolis_hastings, posterior_summary, PriorSpec


def test_map_estimate_matches_mle_with_flat_prior():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=500, random_state=42)
    fitted = map_estimate(data, "gamma")  # priors=None -> plain MLE
    assert fitted.base_shape_params[0] > 0
    assert 1.5 < fitted.base_shape_params[0] < 4.0


def test_map_estimate_with_prior():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=300, random_state=1)
    priors = PriorSpec().add(0, lambda v: -0.5 * ((v - 2.5) / 0.1) ** 2)
    fitted = map_estimate(data, "gamma", priors=priors)
    # strong prior should pull estimate close to 2.5
    assert abs(fitted.base_shape_params[0] - 2.5) < 0.5


def test_metropolis_hastings_runs_and_summarises():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=200, random_state=7)
    result = metropolis_hastings(data, "gamma", n_samples=500, burn_in=200,
                                  proposal_scale=0.1, random_state=0)
    assert result["samples"].shape[0] == 500
    assert 0.0 <= result["acceptance_rate"] <= 1.0
    summary = posterior_summary(result)
    assert "shape1" in summary
    mean, lo, hi = summary["shape1"]
    assert lo <= mean <= hi
