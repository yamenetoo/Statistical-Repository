import numpy as np
import pytest
from unit_transform_dist import UnitTransformMixture, make_unit_transform


def test_mixture_fit_bimodal_data():
    rng = np.random.default_rng(3)
    comp1 = make_unit_transform("gamma", base_shape_params=(2.0,))
    comp2 = make_unit_transform("gamma", base_shape_params=(8.0,), base_scale=0.3)
    data1 = comp1.rvs(size=300, random_state=1)
    data2 = comp2.rvs(size=300, random_state=2)
    data = np.concatenate([data1, data2])
    data = np.clip(data, 1e-4, 1 - 1e-4)

    mix = UnitTransformMixture(["gamma", "gamma"])
    mix.fit(data, max_iter=50, n_init=1, random_state=0)

    assert mix.weights_ is not None
    assert len(mix.weights_) == 2
    assert np.isclose(np.sum(mix.weights_), 1.0, atol=1e-6)

    pdf_vals = mix.pdf(np.array([0.3, 0.5, 0.7]))
    assert np.all(pdf_vals >= 0)

    samples = mix.rvs(size=100, random_state=0)
    assert samples.shape == (100,)
    assert np.all((samples > 0) & (samples < 1))

    aic_val = mix.aic(data)
    bic_val = mix.bic(data)
    assert np.isfinite(aic_val)
    assert np.isfinite(bic_val)


def test_mixture_rejects_out_of_range_data():
    mix = UnitTransformMixture(["gamma", "gamma"])
    with pytest.raises(ValueError):
        mix.fit(np.array([0.5, 1.5, 0.3]))
