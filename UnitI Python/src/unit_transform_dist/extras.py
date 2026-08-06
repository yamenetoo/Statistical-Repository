import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import kstest, cramervonmises
from scipy.optimize import minimize, minimize_scalar, brentq
from scipy.special import chdtri, gamma as gamma_fn
import logging
from typing import Union, Optional, Tuple, Dict, Any

from .core import UnitTransformDist, FittingError, make_unit_transform, fit_unit_transform

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------
# Log‑likelihood, AIC, BIC (with logging)
# -------------------------------------------------------------------

def log_likelihood(dist: UnitTransformDist, data: npt.ArrayLike) -> float:
    data = np.asarray(data)
    return np.sum(dist.logpdf(data))

def aic(dist: UnitTransformDist, data: npt.ArrayLike, num_params: Optional[int] = None) -> float:
    if num_params is None:
        base = dist.base_dist
        num_params = base.numargs + 2
    ll = log_likelihood(dist, data)
    return 2 * num_params - 2 * ll

def bic(dist: UnitTransformDist, data: npt.ArrayLike, num_params: Optional[int] = None) -> float:
    if num_params is None:
        base = dist.base_dist
        num_params = base.numargs + 2
    ll = log_likelihood(dist, data)
    n = len(data)
    return num_params * np.log(n) - 2 * ll

# -------------------------------------------------------------------
# Goodness‑of‑fit tests
# -------------------------------------------------------------------

def ks_test(dist: UnitTransformDist, data: npt.ArrayLike) -> Tuple[float, float]:
    data = np.asarray(data)
    return kstest(data, dist.cdf)

def cramer_von_mises_test(dist: UnitTransformDist, data: npt.ArrayLike):
    data = np.asarray(data)
    u = dist.cdf(data)
    return cramervonmises(u, 'uniform')

def anderson_darling_test(dist: UnitTransformDist, data: npt.ArrayLike) -> Dict[str, Any]:
    data = np.asarray(data)
    u = np.sort(dist.cdf(data))
    n = len(u)
    i = np.arange(1, n + 1)
    log_u = np.log(u)
    log_1_minus_u = np.log(1 - u[::-1])
    A2 = -n - (1 / n) * np.sum((2 * i - 1) * (log_u + log_1_minus_u))
    critical_values = np.array([0.576, 0.656, 0.787, 0.918, 1.092])
    significance_levels = np.array([0.15, 0.10, 0.05, 0.025, 0.01])
    return {
        'statistic': A2,
        'critical_values': critical_values,
        'significance_levels': significance_levels,
    }

# -------------------------------------------------------------------
# Hazard functions
# -------------------------------------------------------------------

def hazard(dist: UnitTransformDist, y: npt.ArrayLike) -> np.ndarray:
    y = np.asarray(y)
    pdf_val = dist.pdf(y)
    sf_val = dist.sf(y)
    with np.errstate(divide='ignore', invalid='ignore'):
        h = pdf_val / sf_val
        h = np.where(sf_val <= 0, np.inf, h)
        h = np.where((y <= 0) | (y >= 1), np.nan, h)
    return h

def cumulative_hazard(dist: UnitTransformDist, y: npt.ArrayLike) -> np.ndarray:
    y = np.asarray(y)
    sf_val = dist.sf(y)
    with np.errstate(divide='ignore'):
        H = -np.log(sf_val)
        H = np.where(sf_val <= 0, np.inf, H)
        H = np.where((y <= 0) | (y >= 1), np.nan, H)
    return H

# -------------------------------------------------------------------
# Q‑Q and P‑P plots
# -------------------------------------------------------------------

def qq_plot(dist: UnitTransformDist, data: npt.ArrayLike, ax=None, show=True, **kwargs):
    data = np.asarray(data)
    if ax is None:
        fig, ax = plt.subplots()
    n = len(data)
    sorted_data = np.sort(data)
    probs = (np.arange(1, n + 1) - 0.5) / n
    theo_quantiles = dist.ppf(probs)
    ax.scatter(theo_quantiles, sorted_data, **kwargs)
    min_val = min(theo_quantiles.min(), sorted_data.min())
    max_val = max(theo_quantiles.max(), sorted_data.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)
    ax.set_xlabel('Theoretical quantiles')
    ax.set_ylabel('Sample quantiles')
    ax.set_title(f'Q‑Q plot for {dist.name}')
    ax.grid(True, alpha=0.3)
    if show:
        plt.show()
    return ax

def pp_plot(dist: UnitTransformDist, data: npt.ArrayLike, ax=None, show=True, **kwargs):
    data = np.asarray(data)
    if ax is None:
        fig, ax = plt.subplots()
    n = len(data)
    sorted_data = np.sort(data)
    emp_cdf = (np.arange(1, n + 1) - 0.5) / n
    theo_cdf = dist.cdf(sorted_data)
    ax.scatter(theo_cdf, emp_cdf, **kwargs)
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
    ax.set_xlabel('Theoretical CDF')
    ax.set_ylabel('Empirical CDF')
    ax.set_title(f'P‑P plot for {dist.name}')
    ax.grid(True, alpha=0.3)
    if show:
        plt.show()
    return ax

# -------------------------------------------------------------------
# Summary report – now using logging and table format
# -------------------------------------------------------------------

def summary(dist: UnitTransformDist, data: npt.ArrayLike, num_params: Optional[int] = None,
            alpha: float = 0.05, verbose: bool = True) -> str:
    data = np.asarray(data)
    n = len(data)
    if num_params is None:
        base = dist.base_dist
        num_params = base.numargs + 2

    ll = log_likelihood(dist, data)
    aic_val = aic(dist, data, num_params)
    bic_val = bic(dist, data, num_params)

    ks_stat, ks_p = ks_test(dist, data)
    try:
        cvm_res = cramer_von_mises_test(dist, data)
        cvm_stat = cvm_res.statistic
        cvm_p = cvm_res.pvalue
    except Exception:
        cvm_stat = cvm_p = np.nan

    ad_res = anderson_darling_test(dist, data)

    rows = []
    rows.append(("Distribution", dist.name))
    rows.append(("Number of observations", n))
    rows.append(("Number of parameters", num_params))
    rows.append(("Log‑likelihood", f"{ll:.4f}"))
    rows.append(("AIC", f"{aic_val:.4f}"))
    rows.append(("BIC", f"{bic_val:.4f}"))
    shape_params = dist.base_shape_params
    if shape_params:
        for i, val in enumerate(shape_params):
            rows.append((f"Shape {i+1}", f"{val:.6f}"))
    rows.append(("Location", f"{dist.base_loc:.6f}"))
    rows.append(("Scale", f"{dist.base_scale:.6f}"))
    rows.append(("KS statistic", f"{ks_stat:.4f}"))
    rows.append(("KS p‑value", f"{ks_p:.4f}"))
    if not np.isnan(cvm_stat):
        rows.append(("CvM statistic", f"{cvm_stat:.4f}"))
        rows.append(("CvM p‑value", f"{cvm_p:.4f}"))
    rows.append(("AD statistic", f"{ad_res['statistic']:.4f}"))
    cv_str = ", ".join([f"{sig:.2f} → {crit:.3f}" for sig, crit in zip(ad_res['significance_levels'], ad_res['critical_values'])])
    rows.append(("AD critical values", cv_str))

    max_key_len = max(len(row[0]) for row in rows) + 2
    max_val_len = max(len(str(row[1])) for row in rows) + 2
    max_key_len = max(max_key_len, 20)
    max_val_len = max(max_val_len, 30)

    separator = "+" + "-" * max_key_len + "+" + "-" * max_val_len + "+"
    header = f"|{'Summary':^{max_key_len}}|{'Value':^{max_val_len}}|"
    lines = [separator, header, separator.replace('-', '=')]
    for key, val in rows:
        lines.append(f"|{str(key):<{max_key_len}}|{str(val):<{max_val_len}}|")
    lines.append(separator)
    report = "\n".join(lines)

    logger.debug(report)
    if verbose:
        print(report)
    return report

# -------------------------------------------------------------------
# Bootstrap confidence intervals (unchanged, but with logging)
# -------------------------------------------------------------------

def bootstrap_ci(dist: UnitTransformDist, data: npt.ArrayLike, n_bootstrap: int = 1000,
                 alpha: float = 0.05, random_state: Optional[Any] = None) -> Dict[str, Tuple[float, float]]:
    data = np.asarray(data)
    n = len(data)
    rng = np.random.default_rng(random_state)
    base_dist = dist.base_dist
    n_shape = base_dist.numargs
    param_names = [f"shape{i+1}" for i in range(n_shape)] + ["loc", "scale"]

    all_params = []
    for i in range(n_bootstrap):
        boot_sample = rng.choice(data, size=n, replace=True)
        try:
            temp = dist.__class__(base_dist=base_dist, name="temp")
            params = temp.fit(boot_sample)
            all_params.append(params)
        except Exception as e:
            logger.warning(f"Bootstrap iteration {i+1} failed: {e}")
            continue

    if not all_params:
        raise RuntimeError("Bootstrap failed for all resamples.")

    all_params = np.array(all_params)
    lower = np.percentile(all_params, 100 * alpha / 2, axis=0)
    upper = np.percentile(all_params, 100 * (1 - alpha / 2), axis=0)

    intervals = {}
    for i, name in enumerate(param_names):
        intervals[name] = (lower[i], upper[i])
    return intervals

# -------------------------------------------------------------------
# Profile likelihood confidence intervals
# -------------------------------------------------------------------

def profile_likelihood_ci(dist: UnitTransformDist, data: npt.ArrayLike,
                          param_index: int, alpha: float = 0.05,
                          bounds: Optional[Tuple[float, float]] = None,
                          max_steps: int = 20) -> Tuple[float, float]:
    """
    Compute confidence interval for a parameter using profile likelihood.
    """
    data = np.asarray(data)
    base_dist = dist.base_dist
    n_shape = base_dist.numargs

    # Fit MLE
    temp = UnitTransformDist(base_dist=base_dist, name="temp")
    mle_params = temp.fit(data)
    shape_params = tuple(mle_params[:n_shape])
    loc = mle_params[-2]
    scale = mle_params[-1]
    fitted = UnitTransformDist(base_dist=base_dist,
                               base_shape_params=shape_params,
                               base_loc=loc, base_scale=scale,
                               name="fitted_temp")
    mle_ll = log_likelihood(fitted, data)
    threshold = 0.5 * chdtri(1, alpha)

    n_params = n_shape + 2
    nuisance_idx = [i for i in range(n_params) if i != param_index]

    # True profile log-likelihood: for a fixed value of the parameter of
    # interest, re-optimise (maximise) over all remaining (nuisance)
    # parameters, rather than merely plugging in their full-MLE values.
    # Fixing nuisance parameters at their unconstrained MLE (as opposed to
    # re-optimising them) understates uncertainty and yields CIs that are
    # too narrow / biased, especially when parameters are correlated
    # (as shape and scale typically are here).
    def profile_ll(p_val):
        if len(nuisance_idx) == 0:
            params = list(mle_params)
            params[param_index] = p_val
        else:
            nuisance_x0 = np.array([mle_params[i] for i in nuisance_idx])

            def neg_ll_nuisance(nuisance_vals):
                params = np.empty(n_params)
                params[param_index] = p_val
                for idx, val in zip(nuisance_idx, nuisance_vals):
                    params[idx] = val
                shp = tuple(params[:n_shape])
                lc, scl = params[-2], params[-1]
                if scl <= 0 or any(s <= 0 for s in shp):
                    return 1e10
                try:
                    td = UnitTransformDist(base_dist=base_dist, base_shape_params=shp,
                                           base_loc=lc, base_scale=scl, name="temp_profile")
                    ll = log_likelihood(td, data)
                    return -ll if np.isfinite(ll) else 1e10
                except Exception:
                    return 1e10

            nuisance_bounds = []
            for idx in nuisance_idx:
                if idx < n_shape or idx == n_params - 1:  # shape or scale: positive
                    nuisance_bounds.append((1e-8, None))
                else:  # loc
                    nuisance_bounds.append((None, None))

            res = minimize(neg_ll_nuisance, nuisance_x0, bounds=nuisance_bounds, method='L-BFGS-B')
            params = np.empty(n_params)
            params[param_index] = p_val
            for idx, val in zip(nuisance_idx, res.x):
                params[idx] = val

        try:
            shape_params = tuple(params[:n_shape])
            loc = params[-2]
            scale = params[-1]
            temp_dist = UnitTransformDist(
                base_dist=base_dist,
                base_shape_params=shape_params,
                base_loc=loc,
                base_scale=scale,
                name="temp"
            )
            return log_likelihood(temp_dist, data)
        except Exception:
            return -np.inf

    # Define function for root‑finding: difference between profile LL and threshold
    def root_func(p):
        return profile_ll(p) - (mle_ll - threshold)

    # Set bounds
    if bounds is None:
        if param_index < n_shape:   # shape parameter (positive)
            bnd = (1e-8, 20.0)      # widened
        elif param_index == n_shape:  # loc
            bnd = (-10.0, 10.0)
        else:                       # scale
            bnd = (1e-8, 20.0)
    else:
        bnd = bounds

    x0 = mle_params[param_index]

    # Find lower bound – search left of MLE
    lower = np.nan
    # Try to find a bracket where root_func changes sign
    left = max(bnd[0], x0 - 0.5)
    # If root_func(left) > 0 (i.e., profile LL still above threshold), move further left
    # until we go below threshold
    step = 0.5
    for _ in range(max_steps):
        if root_func(left) < 0:
            break
        left = max(bnd[0], left - step)
        step *= 1.5
    else:
        # If we couldn't find a left bracket, use a fixed left bound
        left = bnd[0]

    right = x0 + 0.5
    step = 0.5
    for _ in range(max_steps):
        if root_func(right) < 0:
            break
        right = min(bnd[1], right + step)
        step *= 1.5
    else:
        right = bnd[1]

    # Now use bisection or brentq to find roots
    try:
        if left >= bnd[0] and root_func(left) * root_func(x0) <= 0:
            lower = brentq(root_func, left, x0)
        else:
            # fallback: use minimization of absolute value
            res = minimize_scalar(lambda p: abs(root_func(p)),
                                  bounds=(bnd[0], x0),
                                  method='bounded')
            if res.success and res.fun < 1e-4:
                lower = res.x
            else:
                lower = np.nan
    except Exception:
        lower = np.nan

    try:
        if right <= bnd[1] and root_func(x0) * root_func(right) <= 0:
            upper = brentq(root_func, x0, right)
        else:
            res = minimize_scalar(lambda p: abs(root_func(p)),
                                  bounds=(x0, bnd[1]),
                                  method='bounded')
            if res.success and res.fun < 1e-4:
                upper = res.x
            else:
                upper = np.nan
    except Exception:
        upper = np.nan

    return (lower, upper)

# -------------------------------------------------------------------
# Method of Moments (MoM) for common distributions
# -------------------------------------------------------------------

def fit_moments(data: npt.ArrayLike, base_dist_name: str) -> UnitTransformDist:
    """
    Fit using Method of Moments for supported distributions.
    Falls back to MLE if not supported.
    """
    data = np.asarray(data)
    # Transform to X = 1/Y - 1
    x = 1.0 / data - 1.0
    mean = np.mean(x)
    var = np.var(x, ddof=1)  # sample variance

    base_dist = getattr(stats, base_dist_name)
    n_shape = base_dist.numargs

    # We'll support: expon, gamma, lognorm, weibull_min
    if base_dist_name == "expon":
        scale = mean
        shape_params = ()
        loc = 0.0  # exponential is often loc=0
    elif base_dist_name == "gamma":
        if mean > 0 and var > 0:
            shape = mean**2 / var
            scale = var / mean
        else:
            shape = 1.0
            scale = 1.0
        shape_params = (shape,)
        loc = 0.0
    elif base_dist_name == "lognorm":
        # For lognormal: mean and var of log(X) are needed.
        # We approximate using sample log moments.
        logx = np.log(x)
        mu_log = np.mean(logx)
        sigma_log = np.std(logx, ddof=1)
        shape_params = (sigma_log,)
        loc = 0.0
        scale = np.exp(mu_log)  # note: scale = exp(mu)
        # But careful: lognorm in scipy has shape = s, loc, scale.
    elif base_dist_name == "weibull_min":
        # Approximate using mean and variance (requires numerical solving)
        # We'll use a simple approach: set shape=2 (Rayleigh-like) and scale from mean.
        # For a proper MoM, we'd solve two equations.
        shape = 2.0
        # For Weibull, mean = scale * Gamma(1 + 1/shape)
        # => scale = mean / Gamma(1 + 1/shape)
        scale = mean / gamma_fn(1 + 1 / shape)
        shape_params = (shape,)
        loc = 0.0
    else:
        # Fallback to MLE
        logger.warning(f"MoM not implemented for {base_dist_name}, using MLE.")
        return fit_unit_transform(data, base_dist_name)

    # Build the UnitTransformDist with these parameters
    return UnitTransformDist(
        base_dist=base_dist,
        base_shape_params=shape_params,
        base_loc=loc,
        base_scale=scale,
        name=f"unit_{base_dist_name}_mom"
    )

# -------------------------------------------------------------------
# Extended wrapper class (unchanged)
# -------------------------------------------------------------------

class UnitTransformExtras:
    def __init__(self, dist: UnitTransformDist, data: npt.ArrayLike):
        self.dist = dist
        self.data = np.asarray(data)

    def log_likelihood(self):
        return log_likelihood(self.dist, self.data)

    def aic(self, num_params=None):
        return aic(self.dist, self.data, num_params)

    def bic(self, num_params=None):
        return bic(self.dist, self.data, num_params)

    def ks_test(self):
        return ks_test(self.dist, self.data)

    def cramer_von_mises_test(self):
        return cramer_von_mises_test(self.dist, self.data)

    def anderson_darling_test(self):
        return anderson_darling_test(self.dist, self.data)

    def qq_plot(self, **kwargs):
        return qq_plot(self.dist, self.data, **kwargs)

    def pp_plot(self, **kwargs):
        return pp_plot(self.dist, self.data, **kwargs)

    def summary(self, **kwargs):
        return summary(self.dist, self.data, **kwargs)

    def bootstrap_ci(self, **kwargs):
        return bootstrap_ci(self.dist, self.data, **kwargs)

    def profile_ci(self, param_index, **kwargs):
        return profile_likelihood_ci(self.dist, self.data, param_index, **kwargs)

    def hazard(self, y):
        return hazard(self.dist, y)

    def cumulative_hazard(self, y):
        return cumulative_hazard(self.dist, y)

    def __getattr__(self, name):
        if hasattr(self.dist, name):
            return getattr(self.dist, name)
        raise AttributeError(f"'{type(self).__name__}' has no attribute '{name}'")