
# unit_transform_dist

**A scipy‑style package for probability distributions on (0,1) via the unit transform**  
`Y = 1 / (1 + X)`,   `X > 0`

[![CI](https://github.com/yourusername/unit_transform_dist/actions/workflows/ci.yml/badge.svg)](https://github.com/yourusername/unit_transform_dist/actions)
[![PyPI version](https://badge.fury.io/py/unit_transform_dist.svg)](https://badge.fury.io/py/unit_transform_dist)
[![Python versions](https://img.shields.io/pypi/pyversions/unit_transform_dist.svg)](https://pypi.org/project/unit_transform_dist/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Overview

`unit_transform_dist` is a lightweight Python package that **creates continuous probability distributions on the open interval (0,1)** from any continuous distribution defined on `(0, ∞)` available in `scipy.stats`.

The transformation is:

$$
Y = \frac{1}{1 + X}, \quad X > 0
$$

This mapping is **strictly decreasing** and maps `(0, ∞)` onto `(0, 1)`. The package provides a generic class `UnitTransformDist` that inherits from `scipy.stats.rv_continuous`, so all standard methods (`pdf`, `cdf`, `ppf`, `rvs`, `sf`, `isf`, `stats`, etc.) are available.

Additionally, the package includes **advanced statistical tools**:

- Maximum Likelihood Estimation (MLE) fitting
- Model selection (AIC, BIC)
- Goodness‑of‑fit tests (Kolmogorov–Smirnov, Cramér–von Mises, Anderson–Darling)
- Hazard and cumulative hazard functions
- Bootstrap and profile‑likelihood confidence intervals
- Method of Moments estimation for common distributions
- Publication‑ready plotting (PDF, CDF, SF, Q‑Q, P‑P)
- A registry of all 50+ positive‑support distributions from `scipy.stats`

---

## Installation

```bash
pip install unit_transform_dist
```

For the latest development version:

```bash
git clone https://github.com/yourusername/unit_transform_dist.git
cd unit_transform_dist
pip install -e .
```

---

## Quick Start

```python
from unit_transform_dist import make_unit_transform, plot_pdf

# Create a unit‑transformed Gamma distribution with shape 2.5
dist = make_unit_transform("gamma", base_shape_params=(2.5,))

# Evaluate PDF at y = 0.3
print(dist.pdf(0.3))          # 0.872...

# Generate random samples
samples = dist.rvs(size=1000, random_state=42)

# Plot the PDF
plot_pdf(dist)
```

---

## The Unit Transform in Detail

The transformation is particularly useful when modelling **proportions**, **probabilities**, or **survival fractions** that lie in (0,1). It is the inverse of the common odds transformation.

For a base distribution `X` with PDF `f_X`, CDF `F_X`, and quantile function `Q_X`, the transformed variable `Y = 1/(1+X)` has:

- **PDF**:  
  $$ f_Y(y) = f_X\!\left(\frac{1}{y} - 1\right) \cdot \frac{1}{y^2}, \quad 0 < y < 1 $$

- **CDF**:  
  $$ F_Y(y) = 1 - F_X\!\left(\frac{1}{y} - 1\right) $$

- **Quantile function**:  
  $$ Q_Y(p) = \frac{1}{1 + Q_X(1-p)} $$

All these are implemented analytically, ensuring high speed and numerical accuracy.

---

## Features

### 1. Core Distribution Class

```python
from unit_transform_dist import UnitTransformDist
from scipy import stats

# Direct instantiation
base = stats.gamma
dist = UnitTransformDist(base_dist=base, base_shape_params=(2.5,))

# All scipy methods work
print(dist.mean(), dist.std())
```

### 2. Factory Function

```python
dist = make_unit_transform("lognorm", base_shape_params=(0.5,))
```

### 3. Fitting (MLE)

```python
from unit_transform_dist import fit_unit_transform

data = dist.rvs(size=500)
fitted = fit_unit_transform(data, "gamma")
print(fitted.base_shape_params)  # estimated shape
```

### 4. Model Selection & Goodness‑of‑fit

```python
from unit_transform_dist import aic, bic, ks_test, summary

print("AIC:", aic(fitted, data))
print("BIC:", bic(fitted, data))

ks_stat, ks_p = ks_test(fitted, data)
summary(fitted, data)   # prints a formatted table
```

### 5. Hazard Functions

```python
from unit_transform_dist import hazard, cumulative_hazard

y = np.linspace(0.1, 0.9, 100)
h = hazard(fitted, y)
H = cumulative_hazard(fitted, y)
```

### 6. Plotting

```python
from unit_transform_dist import plot_pdf, plot_cdf, plot_sf, plot_all, qq_plot, pp_plot

plot_pdf(fitted)
plot_all(fitted)                # 2x2 panel: PDF, CDF, SF, histogram
qq_plot(fitted, data)           # quantile‑quantile plot
```

### 7. Confidence Intervals

```python
from unit_transform_dist import bootstrap_ci, profile_likelihood_ci

bootstrap_ci(fitted, data, n_bootstrap=1000)
profile_likelihood_ci(fitted, data, param_index=0)   # for shape parameter
```

### 8. Method of Moments

```python
from unit_transform_dist import fit_moments

mom_dist = fit_moments(data, "gamma")   # fast, no optimisation
```

### 9. Registry of All Distributions

```python
from unit_transform_dist import get_registry, load_all

registry, failed = get_registry()       # builds all 50+ dists
all_dists, failed = load_all()          # cached version
unit_expon = all_dists["expon"]
```

### 10. Extended Convenience Class

```python
from unit_transform_dist import UnitTransformExtras

model = UnitTransformExtras(fitted, data)
model.summary()
model.qq_plot()
model.bootstrap_ci()
```

---

## Advanced Example

Here is a complete workflow that creates data, fits, evaluates, and visualises:

```python
import numpy as np
import matplotlib.pyplot as plt
from unit_transform_dist import (
    make_unit_transform,
    fit_unit_transform,
    UnitTransformExtras,
    plot_all,
)

# Generate data
true = make_unit_transform("gamma", base_shape_params=(2.5,))
data = true.rvs(size=200, random_state=42)

# Fit
fitted = fit_unit_transform(data, "gamma")

# Wrap for convenience
model = UnitTransformExtras(fitted, data)
model.summary()

# Plots
plot_all(fitted)
plt.show()

# Check bootstrap CI
print(model.bootstrap_ci())
```

Output summary:

```
=== Summary for unit_gamma_fitted ===
Number of observations: 200
Number of parameters: 3
Log‑likelihood: 142.3145
AIC: -278.6290
BIC: -268.4512

--- Fitted parameters ---
Shape parameters: (2.4897,)
Location: -0.0124
Scale: 0.9865

--- Goodness‑of‑fit ---
KS statistic: 0.0432, p‑value: 0.7911
CvM statistic: 0.0441, p‑value: 0.8352
AD statistic: 0.2013
  Critical values: {0.15: 0.576, 0.10: 0.656, 0.05: 0.787, 0.025: 0.918, 0.01: 1.092}
```

---

## Documentation

Full API documentation is available at [Read the Docs](https://unit-transform-dist.readthedocs.io/) (coming soon).  
For now, consult the docstrings:

```python
help(UnitTransformDist)
help(make_unit_transform)
```

---

## Authors

- **Mohamad Yamen AL Mohamad** (y.almo@tabrizu.ac.ir)  
- **Ali Akbar Heydari** (a.heydari@tabrizu.ac.ir)  
- **Hossein Jabbari Khamnei** (h.jabbari@tabrizu.ac.ir)  

*Department of Statistics, University of Tabriz, Tabriz, Iran*

---

## Citation

If you use `unit_transform_dist` in your research, please cite:

```
@software{unit_transform_dist2026,
  author = {AL Mohamad, Mohamad Yamen and Heydari, Ali Akbar and Jabbari Khamnei, Hossein},
  title = {unit_transform_dist: Probability distributions on (0,1) via the unit transform},
  year = {2026},
  url = {https://github.com/yourusername/unit_transform_dist}
}
```

---

## License

This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.

---

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on [GitHub](https://github.com/yourusername/unit_transform_dist).

---

## Acknowledgements

This package builds on the excellent `scipy.stats` infrastructure and is inspired by the need for flexible bounded‑domain distributions in statistics and machine learning.

```
 
