from typing import Tuple, Union, Optional, Any, Dict
import numpy as np
import numpy.typing as npt
from scipy import stats
from scipy.stats import rv_continuous
from scipy.optimize import minimize
import logging

logger = logging.getLogger(__name__)

class FittingError(Exception):
    pass

class UnitTransformDist(rv_continuous):
    def __init__(self, base_dist: Any, base_shape_params: Tuple[float, ...] = (),
                 base_loc: float = 0.0, base_scale: float = 1.0,
                 name: str = "unit_transform", **kwargs):
        kwargs.setdefault("a", 0.0)
        kwargs.setdefault("b", 1.0)
        kwargs.setdefault("name", name)
        super().__init__(**kwargs)
        self.base_dist = base_dist
        self.base_shape_params = tuple(base_shape_params)
        self.base_loc = base_loc
        self.base_scale = base_scale

    def _base_pdf(self, x: npt.ArrayLike) -> np.ndarray:
        return self.base_dist.pdf(x, *self.base_shape_params,
                                   loc=self.base_loc, scale=self.base_scale)

    def _base_cdf(self, x: npt.ArrayLike) -> np.ndarray:
        return self.base_dist.cdf(x, *self.base_shape_params,
                                   loc=self.base_loc, scale=self.base_scale)

    def _base_ppf(self, q: npt.ArrayLike) -> np.ndarray:
        return self.base_dist.ppf(q, *self.base_shape_params,
                                   loc=self.base_loc, scale=self.base_scale)

    def _pdf(self, y: npt.ArrayLike) -> np.ndarray:
        y = np.asarray(y, dtype=float)
        x = 1.0 / y - 1.0
        fx = self._base_pdf(x)
        jacobian = 1.0 / (y ** 2)
        out = fx * jacobian
        out = np.where((y > 0) & (y < 1), out, 0.0)
        return out

    def _cdf(self, y: npt.ArrayLike) -> np.ndarray:
        y = np.asarray(y, dtype=float)
        x = 1.0 / y - 1.0
        out = 1.0 - self._base_cdf(x)
        return np.clip(out, 0.0, 1.0)

    def _sf(self, y: npt.ArrayLike) -> np.ndarray:
        y = np.asarray(y, dtype=float)
        x = 1.0 / y - 1.0
        return self._base_cdf(x)

    def _ppf(self, q: npt.ArrayLike) -> np.ndarray:
        q = np.asarray(q, dtype=float)
        x_q = self._base_ppf(1.0 - q)
        return 1.0 / (1.0 + x_q)

    def _isf(self, q: npt.ArrayLike) -> np.ndarray:
        q = np.asarray(q, dtype=float)
        x_q = self._base_ppf(q)
        return 1.0 / (1.0 + x_q)

    def _rvs(self, *args, size: Optional[int] = None, random_state: Optional[Any] = None) -> np.ndarray:
        x = self.base_dist.rvs(*self.base_shape_params,
                                loc=self.base_loc, scale=self.base_scale,
                                size=size, random_state=random_state)
        return 1.0 / (1.0 + x)

    def fit(self, data: npt.ArrayLike, **kwargs) -> Tuple[float, ...]:
        data = np.asarray(data)
        if np.any((data <= 0) | (data >= 1)):
            raise ValueError("All data must be in the open interval (0,1).")
        x = 1.0 / data - 1.0
        if hasattr(self.base_dist, 'fit'):
            try:
                return self.base_dist.fit(x, **kwargs)
            except Exception as e:
                raise FittingError(f"Base distribution fit failed: {e}")
        # Generic optimizer
        num_params = self.base_dist.numargs + 2
        x_std = np.std(x) or 1.0
        x0 = np.ones(num_params)
        x0[-2] = 0.0
        x0[-1] = x_std
        bounds = [(1e-8, None)] * self.base_dist.numargs + [(None, None), (1e-8, None)]

        def neg_log_likelihood(params):
            shape_params = tuple(params[:-2])
            loc = params[-2]
            scale = params[-1]
            if scale <= 0:
                return 1e10
            try:
                logpdf = self.base_dist.logpdf(x, *shape_params, loc=loc, scale=scale)
                return -np.sum(logpdf)
            except Exception:
                return 1e10

        result = minimize(neg_log_likelihood, x0, bounds=bounds, method='L-BFGS-B')
        if not result.success:
            raise FittingError(f"Optimization failed: {result.message}")
        return tuple(result.x)

def make_unit_transform(base_dist_name: str, base_shape_params: Tuple[float, ...] = (),
                        base_loc: float = 0.0, base_scale: float = 1.0) -> UnitTransformDist:
    base_dist = getattr(stats, base_dist_name)
    name = f"unit_{base_dist_name}"
    return UnitTransformDist(base_dist=base_dist,
                              base_shape_params=base_shape_params,
                              base_loc=base_loc,
                              base_scale=base_scale,
                              name=name)

def fit_unit_transform(data: npt.ArrayLike, base_dist_name: str, **kwargs) -> UnitTransformDist:
    base_dist = getattr(stats, base_dist_name)
    temp_dist = UnitTransformDist(base_dist=base_dist, name=f"unit_{base_dist_name}")
    params = temp_dist.fit(data, **kwargs)
    n_shape = base_dist.numargs
    shape_params = params[:n_shape]
    loc = params[-2] if len(params) >= 2 else 0.0
    scale = params[-1] if len(params) >= 1 else 1.0
    return UnitTransformDist(
        base_dist=base_dist,
        base_shape_params=shape_params,
        base_loc=loc,
        base_scale=scale,
        name=f"unit_{base_dist_name}_fitted"
    )