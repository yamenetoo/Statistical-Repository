import numpy as np
import pytest
from unit_transform_dist import UnitTransformRegression, UnitTransformDist
from scipy import stats


def test_regression_fits_and_predicts():
    rng = np.random.default_rng(0)
    n = 300
    x1 = rng.normal(size=n)
    true_beta = np.array([0.2, 0.5])  # intercept, slope on log-scale link
    eta = true_beta[0] + true_beta[1] * x1
    scale = np.exp(eta)
    shape = (2.0,)

    y = np.array([
        UnitTransformDist(base_dist=stats.gamma, base_shape_params=shape,
                           base_loc=0.0, base_scale=s).rvs(size=1, random_state=rng)[0]
        for s in scale
    ])

    model = UnitTransformRegression("gamma", link="log", fit_shape=True)
    model.fit(x1.reshape(-1, 1), y)

    assert model.coef_ is not None
    assert len(model.coef_) == 2
    preds = model.predict_scale(x1.reshape(-1, 1))
    assert preds.shape == (n,)
    assert np.all(preds > 0)

    means = model.predict_mean(x1.reshape(-1, 1))
    assert means.shape == (n,)
    assert np.all((means > 0) & (means < 1))

    ll = model.loglik()
    assert np.isfinite(ll)


def test_regression_requires_unit_interval():
    model = UnitTransformRegression("gamma")
    X = np.random.default_rng(0).normal(size=(10, 1))
    y = np.ones(10) * 1.5  # invalid: outside (0,1)
    with pytest.raises(ValueError):
        model.fit(X, y)
