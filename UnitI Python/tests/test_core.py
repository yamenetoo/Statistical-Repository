import pytest
import numpy as np
from unit_transform_dist import make_unit_transform, fit_unit_transform, UnitTransformDist

def test_pdf_integrates_to_one():
    dist = make_unit_transform("expon")
    y = np.linspace(0.001, 0.999, 1000)
    pdf = dist.pdf(y)
    integral = np.trapezoid(pdf, y)
    assert np.isclose(integral, 1.0, rtol=1e-2)

def test_cdf_monotonic():
    dist = make_unit_transform("gamma", base_shape_params=(2.0,))
    y = np.linspace(0.01, 0.99, 100)
    cdf = dist.cdf(y)
    assert np.all(np.diff(cdf) >= 0)

def test_fit_recovers_parameters():
    true = make_unit_transform("gamma", base_shape_params=(2.5,))
    data = true.rvs(size=1000, random_state=42)
    fitted = fit_unit_transform(data, "gamma")
    # Check shape is close
    assert np.isclose(fitted.base_shape_params[0], 2.5, rtol=0.2)