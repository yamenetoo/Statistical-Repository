from .core import UnitTransformDist, make_unit_transform, fit_unit_transform, FittingError
from .generate_all import get_positive_distributions, build_all_unit_transforms, load_all
from .registry import get_registry, register_custom, get_all_distributions
from .plot import plot_pdf, plot_cdf, plot_sf, plot_all
from .extras import (
    log_likelihood, aic, bic, ks_test, cramer_von_mises_test,
    anderson_darling_test, hazard, cumulative_hazard,
    qq_plot, pp_plot, summary, bootstrap_ci,
    profile_likelihood_ci, fit_moments, UnitTransformExtras
)
from .utils import setup_logging
from .lmoments import sample_lmoments, lmom_ratios, fit_lmoments
from .regression import UnitTransformRegression
from .mixture import UnitTransformMixture
from .bayes import map_estimate, metropolis_hastings, posterior_summary, PriorSpec
from .censored import fit_censored
from .diagnostics import (
    pearson_residuals, cox_snell_residuals, quantile_residuals, cooks_distance
)
from .model_selection import auto_select_distribution

__all__ = [
    # core
    "UnitTransformDist", "make_unit_transform", "fit_unit_transform", "FittingError",
    # registry / generation
    "get_positive_distributions", "build_all_unit_transforms", "load_all",
    "get_registry", "register_custom", "get_all_distributions",
    # plotting
    "plot_pdf", "plot_cdf", "plot_sf", "plot_all",
    # goodness-of-fit / summary
    "log_likelihood", "aic", "bic", "ks_test", "cramer_von_mises_test",
    "anderson_darling_test", "hazard", "cumulative_hazard",
    "qq_plot", "pp_plot", "summary", "bootstrap_ci",
    "profile_likelihood_ci", "fit_moments", "UnitTransformExtras",
    "setup_logging",
    # L-moments
    "sample_lmoments", "lmom_ratios", "fit_lmoments",
    # regression
    "UnitTransformRegression",
    # mixtures
    "UnitTransformMixture",
    # Bayesian
    "map_estimate", "metropolis_hastings", "posterior_summary", "PriorSpec",
    # censoring / truncation
    "fit_censored",
    # diagnostics
    "pearson_residuals", "cox_snell_residuals", "quantile_residuals", "cooks_distance",
    # model selection
    "auto_select_distribution",
]

__version__ = "0.3.0"
