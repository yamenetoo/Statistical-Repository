import numpy as np
import pytest
from unit_transform_dist import make_unit_transform, auto_select_distribution


def test_auto_select_recovers_true_family_in_top_results():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=500, random_state=42)
    results = auto_select_distribution(data, candidates=["gamma", "expon", "lognorm", "weibull_min"],
                                        criterion="aic")
    assert len(results) > 0
    top_names = [r["name"] for r in results[:2]]
    assert "gamma" in top_names


def test_auto_select_top_n_and_criteria():
    true = make_unit_transform("expon")
    data = true.rvs(size=300, random_state=1)
    results = auto_select_distribution(data, candidates=["gamma", "expon", "lognorm"],
                                        criterion="bic", top_n=2)
    assert len(results) == 2
    # results sorted ascending by bic
    assert results[0]["bic"] <= results[1]["bic"]


def test_auto_select_invalid_criterion_raises():
    true = make_unit_transform("expon")
    data = true.rvs(size=50, random_state=1)
    with pytest.raises(ValueError):
        auto_select_distribution(data, candidates=["expon"], criterion="not_a_criterion")
