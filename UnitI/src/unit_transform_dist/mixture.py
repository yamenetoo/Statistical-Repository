"""
mixture.py
==========

Finite mixtures of UnitTransformDist components, fitted by the EM
algorithm. Useful for bimodal or multimodal data on (0,1).
"""

from typing import List, Optional, Tuple
import numpy as np
import numpy.typing as npt
from scipy import stats

from .core import UnitTransformDist, FittingError

__all__ = ["UnitTransformMixture"]


class UnitTransformMixture:
    """
    K-component mixture of UnitTransformDist distributions.

    Parameters
    ----------
    base_dist_names : list of str
        Base distribution name for each component (length K). May repeat
        the same name for several components (e.g. two Gammas).
    """

    def __init__(self, base_dist_names: List[str]):
        self.base_dist_names = list(base_dist_names)
        self.k = len(base_dist_names)
        self.base_dists = [getattr(stats, n) for n in base_dist_names]
        self.weights_: Optional[np.ndarray] = None
        self.components_: Optional[List[UnitTransformDist]] = None
        self.n_iter_: int = 0
        self.loglik_history_: List[float] = []

    def _component_pdf_matrix(self, data: np.ndarray) -> np.ndarray:
        mat = np.zeros((len(data), self.k))
        for j, comp in enumerate(self.components_):
            mat[:, j] = comp.pdf(data)
        return mat

    def fit(self, data: npt.ArrayLike, max_iter: int = 200, tol: float = 1e-6,
            random_state: Optional[int] = None, n_init: int = 3) -> "UnitTransformMixture":
        data = np.asarray(data, dtype=float)
        if np.any((data <= 0) | (data >= 1)):
            raise ValueError("All data must be in the open interval (0,1).")

        rng = np.random.default_rng(random_state)
        best_ll = -np.inf
        best_state = None

        for init in range(n_init):
            weights = np.full(self.k, 1.0 / self.k)
            # random component assignment for initial parameter guesses
            components = []
            splits = np.array_split(rng.permutation(data), self.k)
            for j, (base_dist, chunk) in enumerate(zip(self.base_dists, splits)):
                chunk = chunk if len(chunk) >= max(3, base_dist.numargs + 2) else data
                try:
                    temp = UnitTransformDist(base_dist=base_dist, name=f"comp{j}")
                    params = temp.fit(chunk)
                    n_shape = base_dist.numargs
                    shape_params = tuple(params[:n_shape])
                    loc, scale = params[-2], params[-1]
                except Exception:
                    shape_params = tuple(np.ones(base_dist.numargs))
                    loc, scale = 0.0, 1.0
                components.append(UnitTransformDist(base_dist=base_dist, base_shape_params=shape_params,
                                                      base_loc=loc, base_scale=scale, name=f"comp{j}"))

            ll_hist = []
            prev_ll = -np.inf
            for it in range(max_iter):
                self.components_ = components
                pdf_mat = self._component_pdf_matrix(data)
                weighted = pdf_mat * weights[np.newaxis, :]
                row_sums = weighted.sum(axis=1)
                row_sums = np.where(row_sums <= 0, 1e-300, row_sums)
                resp = weighted / row_sums[:, np.newaxis]

                ll = np.sum(np.log(row_sums))
                ll_hist.append(ll)
                if abs(ll - prev_ll) < tol:
                    break
                prev_ll = ll

                # M-step
                nk = resp.sum(axis=0)
                weights = nk / len(data)
                new_components = []
                for j, base_dist in enumerate(self.base_dists):
                    w = resp[:, j]
                    if nk[j] < max(3, base_dist.numargs + 2):
                        new_components.append(components[j])
                        continue
                    try:
                        new_components.append(
                            self._weighted_fit(base_dist, data, w, name=f"comp{j}")
                        )
                    except Exception:
                        new_components.append(components[j])
                components = new_components

            if ll_hist and ll_hist[-1] > best_ll:
                best_ll = ll_hist[-1]
                best_state = (weights.copy(), components, ll_hist, it + 1)

        self.weights_, self.components_, self.loglik_history_, self.n_iter_ = best_state
        return self

    @staticmethod
    def _weighted_fit(base_dist, data: np.ndarray, weights: np.ndarray, name: str) -> UnitTransformDist:
        """Weighted MLE fit of a single component via optimisation of weighted log-likelihood."""
        from scipy.optimize import minimize
        n_shape = base_dist.numargs
        x = 1.0 / data - 1.0
        x0 = np.ones(n_shape + 2)
        x0[-2] = 0.0
        x0[-1] = max(np.average(x, weights=weights), 1e-3)

        def neg_ll(params):
            shape_params = tuple(params[:n_shape])
            loc, scale = params[-2], params[-1]
            if scale <= 0:
                return 1e10
            try:
                logpdf = base_dist.logpdf(x, *shape_params, loc=loc, scale=scale)
                val = -np.sum(weights * logpdf)
                return val if np.isfinite(val) else 1e10
            except Exception:
                return 1e10

        # Bound loc so the fitted base-distribution support does not drift
        # into negative territory (which would push Y = 1/(1+X) outside (0,1)).
        loc_lower = -0.999 * np.min(x)
        bounds = [(1e-8, None)] * n_shape + [(loc_lower, None), (1e-8, None)]
        result = minimize(neg_ll, x0, bounds=bounds, method='L-BFGS-B')
        params = result.x
        shape_params = tuple(params[:n_shape])
        loc, scale = params[-2], params[-1]
        return UnitTransformDist(base_dist=base_dist, base_shape_params=shape_params,
                                  base_loc=loc, base_scale=scale, name=name)

    def pdf(self, y: npt.ArrayLike) -> np.ndarray:
        y = np.asarray(y, dtype=float)
        mat = self._component_pdf_matrix(y)
        return mat @ self.weights_

    def cdf(self, y: npt.ArrayLike) -> np.ndarray:
        y = np.asarray(y, dtype=float)
        mat = np.column_stack([comp.cdf(y) for comp in self.components_])
        return mat @ self.weights_

    def rvs(self, size: int = 1, random_state: Optional[int] = None) -> np.ndarray:
        rng = np.random.default_rng(random_state)
        comp_idx = rng.choice(self.k, size=size, p=self.weights_)
        out = np.empty(size)
        for j in range(self.k):
            mask = comp_idx == j
            n_j = mask.sum()
            if n_j:
                out[mask] = self.components_[j].rvs(size=n_j, random_state=rng)
        # Defensive clip: guards against rare floating-point/degenerate-fit
        # edge cases where a component's support drifts slightly outside (0,1).
        eps = 1e-10
        return np.clip(out, eps, 1 - eps)

    def aic(self, data: npt.ArrayLike) -> float:
        n_params = sum(c.base_dist.numargs + 2 for c in self.components_) + (self.k - 1)
        ll = np.sum(np.log(self.pdf(data)))
        return 2 * n_params - 2 * ll

    def bic(self, data: npt.ArrayLike) -> float:
        n_params = sum(c.base_dist.numargs + 2 for c in self.components_) + (self.k - 1)
        ll = np.sum(np.log(self.pdf(data)))
        return n_params * np.log(len(data)) - 2 * ll
