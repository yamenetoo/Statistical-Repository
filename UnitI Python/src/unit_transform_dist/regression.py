"""
regression.py
=============

Regression framework for bounded (0,1) responses using UnitTransformDist,
in the spirit of beta regression (Ferrari & Cribari-Neto, 2004). One (or
more) of the base distribution's parameters is linked to a linear predictor
through a monotone link function, and all regression coefficients are
estimated jointly by maximum likelihood.

Currently the scale parameter of the base distribution is modelled as
scale_i = g^{-1}(x_i' beta), with g = log by default (ensuring positivity).
Shape parameters and location can optionally be estimated jointly (assumed
constant across observations).
"""

from typing import Optional, Tuple, Callable
import numpy as np
import numpy.typing as npt
from scipy import stats
from scipy.optimize import minimize

from .core import UnitTransformDist, FittingError

__all__ = ["UnitTransformRegression"]


def _log_link(mu):
    """mu -> eta (linear predictor scale)."""
    return np.log(mu)


def _log_link_inv(eta):
    """eta -> mu (parameter scale, always positive)."""
    return np.exp(eta)


def _identity_link(x):
    return x


_LINKS = {
    # name: (link mu->eta, inverse-link eta->mu)
    "log": (_log_link, _log_link_inv),
    "identity": (_identity_link, _identity_link),
}


class UnitTransformRegression:
    """
    Regression model for (0,1)-valued responses via the unit transform.

    y_i ~ UnitTransformDist(base_dist, shape_params, loc=0, scale=scale_i)
    scale_i = link^{-1}(X_i @ beta)

    Parameters
    ----------
    base_dist_name : str
        Name of a scipy.stats positive-support continuous distribution.
    link : str, default "log"
        Link function for the scale parameter ("log" or "identity").
    fit_shape : bool, default True
        If True, shape parameters of the base distribution are estimated
        jointly (assumed constant across observations). If False, they must
        be supplied via `fixed_shape_params`.
    fixed_shape_params : tuple, optional
        Fixed shape parameters to use when fit_shape=False.

    Examples
    --------
    >>> model = UnitTransformRegression("gamma")
    >>> model.fit(X, y)
    >>> model.predict(X_new)
    """

    def __init__(self, base_dist_name: str, link: str = "log",
                 fit_shape: bool = True, fixed_shape_params: Optional[Tuple[float, ...]] = None):
        if link not in _LINKS:
            raise ValueError(f"Unknown link '{link}'. Choose from {list(_LINKS)}.")
        self.base_dist_name = base_dist_name
        self.base_dist = getattr(stats, base_dist_name)
        self.link_name = link
        self.link_fn, self.link_inv_fn = _LINKS[link]
        self.fit_shape = fit_shape
        self.fixed_shape_params = fixed_shape_params
        self.n_shape = self.base_dist.numargs

        self.coef_: Optional[np.ndarray] = None
        self.shape_params_: Optional[Tuple[float, ...]] = None
        self.n_features_: Optional[int] = None
        self.result_ = None

    def _unpack(self, params: np.ndarray, n_features: int):
        beta = params[:n_features]
        if self.fit_shape:
            shape_params = tuple(params[n_features:n_features + self.n_shape])
        else:
            shape_params = self.fixed_shape_params or ()
        return beta, shape_params

    def _neg_log_likelihood(self, params: np.ndarray, X: np.ndarray, y: np.ndarray) -> float:
        n_features = X.shape[1]
        beta, shape_params = self._unpack(params, n_features)
        eta = X @ beta
        scale = self.link_inv_fn(eta)
        if np.any(scale <= 0) or not np.all(np.isfinite(scale)):
            return 1e10
        try:
            dist = UnitTransformDist(base_dist=self.base_dist, base_shape_params=shape_params,
                                      base_loc=0.0, base_scale=1.0, name="temp_reg")
            # scale enters per-observation: evaluate logpdf with varying scale
            # by rescaling x = 1/y - 1 manually for speed/robustness.
            x = 1.0 / y - 1.0
            logpdf = self.base_dist.logpdf(x, *shape_params, loc=0.0, scale=scale)
            jac = -2.0 * np.log(y)  # log(1/y^2)
            ll = np.sum(logpdf + jac)
            if not np.isfinite(ll):
                return 1e10
            return -ll
        except Exception:
            return 1e10

    def fit(self, X: npt.ArrayLike, y: npt.ArrayLike, add_intercept: bool = True,
            x0: Optional[npt.ArrayLike] = None, maxiter: int = 500):
        X = np.atleast_2d(np.asarray(X, dtype=float))
        if X.shape[0] == 1 and X.shape[1] != 1:
            X = X.T
        y = np.asarray(y, dtype=float)
        if np.any((y <= 0) | (y >= 1)):
            raise ValueError("All responses must lie in the open interval (0,1).")
        if add_intercept:
            X = np.column_stack([np.ones(X.shape[0]), X])
        n_features = X.shape[1]
        self.n_features_ = n_features

        n_params = n_features + (self.n_shape if self.fit_shape else 0)
        if x0 is None:
            x0 = np.zeros(n_params)
            # reasonable starting scale via mean of x = 1/y - 1
            x_transformed = 1.0 / y - 1.0
            start_scale = max(np.mean(x_transformed), 1e-3)
            x0[0] = self.link_fn(start_scale)
            if self.fit_shape:
                x0[n_features:] = 1.0

        bounds = [(None, None)] * n_features
        if self.fit_shape:
            bounds += [(1e-6, None)] * self.n_shape

        result = minimize(self._neg_log_likelihood, x0, args=(X, y),
                           method='L-BFGS-B', bounds=bounds,
                           options={"maxiter": maxiter})
        if not result.success:
            raise FittingError(f"Regression optimisation failed: {result.message}")

        self.coef_, self.shape_params_ = self._unpack(result.x, n_features)
        self.result_ = result
        self._add_intercept = add_intercept
        return self

    def predict_scale(self, X: npt.ArrayLike, add_intercept: Optional[bool] = None) -> np.ndarray:
        if self.coef_ is None:
            raise FittingError("Model has not been fit yet.")
        X = np.atleast_2d(np.asarray(X, dtype=float))
        if X.shape[0] == 1 and X.shape[1] != self.n_features_ - int(self._add_intercept):
            X = X.T
        add_intercept = self._add_intercept if add_intercept is None else add_intercept
        if add_intercept:
            X = np.column_stack([np.ones(X.shape[0]), X])
        eta = X @ self.coef_
        return self.link_inv_fn(eta)

    def predict_mean(self, X: npt.ArrayLike) -> np.ndarray:
        """Predict E[Y | X] via the fitted distribution's mean for each row."""
        scales = self.predict_scale(X)
        means = np.array([
            UnitTransformDist(base_dist=self.base_dist, base_shape_params=self.shape_params_,
                               base_loc=0.0, base_scale=s, name="temp").mean()
            for s in scales
        ])
        return means

    def loglik(self) -> float:
        if self.result_ is None:
            raise FittingError("Model has not been fit yet.")
        return -self.result_.fun

    def summary(self) -> str:
        lines = [f"UnitTransformRegression(base_dist='{self.base_dist_name}', link='{self.link_name}')",
                 f"  Coefficients (scale link): {np.round(self.coef_, 6)}"]
        if self.fit_shape:
            lines.append(f"  Shape parameters: {np.round(self.shape_params_, 6)}")
        lines.append(f"  Log-likelihood: {self.loglik():.4f}")
        text = "\n".join(lines)
        print(text)
        return text
