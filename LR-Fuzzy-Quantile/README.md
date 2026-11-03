# LR-Fuzzy-Quantile

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yamenetoo/Statistical-Repository/blob/main/LR-Fuzzy-Quantile/lr_fuzzy_simulation.ipynb)

Python code accompanying the article

> **Quantile representation of LR random fuzzy numbers**
> Mohamad Yamen Al-Mohamad, Abbas Parchami, Ali Akbar Hydari
> (a contribution to Parchami, Grzegorzewski & Romaniuk (2024), *Statistical simulations with LR random fuzzy numbers*, Statistical Papers 65, 3585-3609).

The repository contains (i) the simulation study of the article, (ii) the inference experiments of Section 4, and (iii) an independent Python implementation of the base paper, so that the two generators can be compared side by side.

---

## 1. What is the contribution?

An **LR random fuzzy number** `X = (C, S^l, S^r)_LR` is built from a random core `C` and two random spreads `S^l, S^r`;
each arm is `1 - F_S(t)/F_S(s)`, the distribution function of the spread truncated to its realised value `s`.
Parchami et al. generate it either as a limit of `k`-knot polygons (Algorithm 1) or from the truncated c.d.f. (Algorithm 2).

The article shows that, when the spread laws have continuous strictly increasing c.d.f.s `F_l, F_r` (quantile functions `Q_l, Q_r`), the limit has the exact **quantile representation**

```
X_alpha = [ C - Q_l((1-alpha) V_l) ,  C + Q_r((1-alpha) V_r) ],     V_j = F_j(S^j) ~ U(0,1),  (C, V_l, V_r) independent.
```

Consequences studied (and simulated here):

| Result | Statement | Simulated in |
|---|---|---|
| Generator (Algorithm 2') | three independent variates and two quantile evaluations per fuzzy number | `lr_fuzzy_simulation.py` (equivalence, timing) |
| Aumann expectation | `E X_alpha = [E C - m_l(alpha), E C + m_r(alpha)]`, `m(alpha) = (1/(1-alpha)) int_0^{1-alpha} Q(u) du` (closed forms for Exp, U, Beta(a,1), Weibull) | `lr_fuzzy_simulation.py` (Aumann sheet) |
| Non-asymptotic rate | `P( sup|X_k - X| <= sqrt(ln(4/delta)/(2k)) + 2/(k+1) ) >= 1 - delta` | `lr_fuzzy_simulation.py` (Convergence sheet) |
| Dependence | Gaussian copula on `(U_C, V_l, V_r)`; mean width unchanged, variance of the width increases with `rho`; fuzzy intervals | `lr_fuzzy_simulation.py` (Dependence sheet) |
| Inference | exact recovery of the spread law at rate `1/n`, Kolmogorov-type goodness-of-fit test, Fisher information of the knots | `lr_fuzzy_inference.py` |
| Horizontal rate | `P( d_inf(X_k, X) <= max_j omega_{Q_j,V_j}(V_j eps_k(delta)) ) >= 1 - delta` (`d_inf` = largest Hausdorff distance between alpha-cuts, `omega` = modulus of continuity of the quantile function); exact order `k^(-1/2)` when `Q` is Lipschitz | `lr_fuzzy_horizontal.py` |
| Characterisation of the LR class | a spread-determined random fuzzy number is LR iff its arms satisfy multiplicative consistency `G_s(t) = G_s(s') G_{s'}(t)` and the shape function so identified is the law of the spread | `lr_fuzzy_horizontal.py` |

---

## 2. Repository layout

```
LR-Fuzzy-Quantile/
├── README.md                       this file
├── requirements.txt                numpy, scipy, pandas, openpyxl
├── lr_fuzzy_simulation.py          main simulation study (writes an Excel workbook)
├── lr_fuzzy_simulation.ipynb       the same study as a Google Colab notebook
├── lr_fuzzy_inference.py           inference experiments of Section 4 (writes a JSON file)
├── lr_fuzzy_horizontal.py          horizontal-rate and characterisation experiments (writes a JSON file)
├── results/horizontal_results.json output of the full run of lr_fuzzy_horizontal.py (seed 2026)
└── base_paper/
    ├── README.md
    ├── parchami2024_base_paper.py  Algorithms 1-3 and Examples 1-3 of the base paper
    └── figures/                    fig1_example1.png, fig3_example2.png, fig4_example3.png
```

---

## 3. Quick start

### Google Colab (no installation)
Click the badge above, then **Runtime -> Run all**. The first cell installs the packages and the last cell downloads the Excel file.

### Local machine
```bash
git clone https://github.com/yamenetoo/Statistical-Repository.git
cd Statistical-Repository/LR-Fuzzy-Quantile
pip install -r requirements.txt            # on some systems: pip install --break-system-packages -r requirements.txt

python lr_fuzzy_simulation.py --quick      # fast test (small numbers of replications)
python lr_fuzzy_simulation.py              # default: 500 replications, seed 2026
python lr_fuzzy_simulation.py --reps 2000 --seed 7 --out my_results.xlsx

python lr_fuzzy_inference.py --quick       # inference experiments, test run
python lr_fuzzy_inference.py               # full run (several minutes)
python lr_fuzzy_inference.py --only gof --out gof.json

python lr_fuzzy_horizontal.py --quick     # horizontal rate + characterisation, test run
python lr_fuzzy_horizontal.py             # full run (about 3-4 minutes)
python lr_fuzzy_horizontal.py --out h.json --datadir figdata   # also writes .dat files for pgfplots

cd base_paper && python parchami2024_base_paper.py --example all --outdir figures
```
Requirements: Python 3.9+ ; numpy, scipy, pandas, openpyxl ; matplotlib for the base-paper figures.

---

## 4. The experiments

### 4.1 `lr_fuzzy_simulation.py` -> one Excel sheet per experiment

| Sheet | Experiment | What to look for |
|---|---|---|
| `Config` | seed, replications, `k` values, `delta`, scenarios | records the settings of the run |
| `Equivalence` | membership of the paper's Algorithm 2 evaluated at the endpoints given by the quantile formula | must equal `alpha`; maximal absolute error is at machine precision (about 1e-15 in the article) |
| `Convergence` | sup-distance between the `k`-knot polygon (Algorithm 1) and the LR limit, `k = 2 ... 1000` | mean error decays like `1/sqrt(k)`; `Coverage_of_bound` is the empirical frequency with which the error is below the bound of the article |
| `Aumann` | numerical integral, closed form and Monte-Carlo mean of the expected `alpha`-cut offsets | differences within Monte-Carlo error (`MC_std_error`) |
| `Dependence` | Gaussian copula between `V_l` and `V_r` (and between `C` and the spreads) | `Mean_width` constant in `rho`; `Var_width` increasing in `rho` |
| `Timing` | Algorithm 1 vs Algorithm 2 vs Algorithm 2' (101 `alpha` levels, 1000 fuzzy numbers) | `Speedup_vs_Alg1`; the article reports about 3x for Algorithm 2' |
| `Fuzzy_sample_*` | a simulated fuzzy random sample (long and wide format) | ready for plotting / further use |
| `Fuzzy_interval_*` | the same with a random core width `W` (fuzzy intervals) | trapezoid-like numbers |

Scenarios (edit `SCENARIOS` at the top of the script to try other laws):

| | Core `C` | Left spread `S^l` | Right spread `S^r` |
|---|---|---|---|
| A | N(0,1) | Exp(3) | Beta(3,1) |
| B | N(0,1) | Chi2(8) | Chi2(8) |
| C | N(0,1) | Gamma(2,1) | U(0,3) |

Command-line options: `--reps` (replications per cell, default 500), `--seed` (default 2026), `--out` (output file), `--quick`.

### 4.2 `lr_fuzzy_inference.py` -> `inference_results.json`

Exp(3) spreads, noise-free observed memberships.

| Block | Result | Quantity |
|---|---|---|
| `recover` | Proposition 4.1 | error of the spread c.d.f. estimated from the largest-spread fuzzy number; `n * error` stays near the limit `1 + exp(-2)` (about 1.14), against `sqrt(n) * error` for the crisp empirical c.d.f. |
| `compat` | Lemma 4.3 | the rescaled-truncation identity holds for LR (and triangular) arms: necessary, not sufficient, for the null hypothesis |
| `gof` | Proposition 4.2 | rejection rates of the Kolmogorov-type test (nominal 5 %) under H0, triangular arms and a misspecified spread law |
| `mle` | Proposition 4.4 | standard deviation of the MLE of `lambda` using `k` knots per arm against the Fisher-information prediction `1/sqrt(n I_k)`, and the gain `I_k / I_0` |

Options: `--seed`, `--quick`, `--only {recover,compat,gof,mle}`, `--out`.

### 4.3 `lr_fuzzy_horizontal.py` -> `horizontal_results.json`

| Block | Result | Quantity |
|---|---|---|
| `horizontal` | Theorem "horizontal bound" | for scenarios A-C and `k = 4 ... 1000`: median / 95% quantile of `d_inf(X_k, X)`, the bound with `delta = 0.05` (modulus of continuity evaluated on a grid), its empirical coverage and the ratio bound/error |
| `fitted_slopes` | Remark on Hoelder exponents | observed log-log order `k^(-slope)` of the median error (about 0.46, 0.34, 0.30 for A, C, B in the stored run) |
| `sharp` | Proposition "exact order" | uniform spread, `V = 1`: `sqrt(k) * E d_inf` tends to `sqrt(pi/2) ln 2 = 0.8687` |
| `characterisation` | Theorem "characterisation" | largest violation of the consistency condition (CE) for LR (exponential) arms, triangular arms and a spread-determined non-LR family, and the error of recovering the shape function from the arms alone |

Options: `--seed`, `--quick`, `--out`, `--datadir`.

### 4.4 `base_paper/parchami2024_base_paper.py`

| Paper | Function |
|---|---|
| Algorithm 1 (k-knot polygon) | `knot_fuzzy_number`, `knot_from_realisation` |
| Algorithm 2 (LR random fuzzy number) | `lr_fuzzy_number`, class `LRFuzzyNumber` (shapes `L`, `R`, eqs. (14)-(15)) |
| Algorithm 3 (fuzzy random sample) | `fuzzy_random_sample` |
| Example 1 / Fig. 1 | `example1` - reproduces the knots (1.689, 0.333), (1.700, 0.667), (1.873, 0.333), (1.769, 0.667) |
| Example 2 / Fig. 3 | `example2` - X_k, Y_k, Z_k for k = 2, 4, 8, 16, 50, 200 and the LR limits |
| Example 3 / Fig. 4 | `example3` - checks `L(x) = 1 - x`, `R(x) = 1 - x^3` |

Conventions: right-arm knots use `c + r_(k-i+1)` as in Example 1 of the paper (the printed eq. (7) has `r_(i)`); `Gamma(1,6)` of Example 2 is read as shape 1, rate 6.

---

## 5. Using the generator in your own code

```python
import numpy as np
from scipy import stats

rng = np.random.default_rng(1)
Q_l = stats.expon(scale=1/3).ppf          # quantile function of S^l
Q_r = stats.beta(3, 1).ppf                # quantile function of S^r
alpha = np.linspace(0, 1, 101)

# one LR random fuzzy number (Algorithm 2'): three uniforms, two quantile evaluations per level
C, Vl, Vr = stats.norm(0, 1).rvs(random_state=rng), rng.random(), rng.random()
left  = C - Q_l((1 - alpha) * Vl)         # left end of the alpha-cuts
right = C + Q_r((1 - alpha) * Vr)         # right end of the alpha-cuts
```
The library functions `quantile_offset`, `exp_sample` (fuzzy random samples, optional core width `W` for fuzzy intervals) and `gaussian_copula_uniforms` in `lr_fuzzy_simulation.py` do the same in vectorised form.

---

## 6. Reproducibility

* Default seeds: 2026 (simulation and inference scripts), 2024 (base-paper figures); results are deterministic for a given NumPy version.
* The numbers in the article were produced with Python 3.13, NumPy 2.5.3 and SciPy 1.18.1 on a virtual machine with two Intel Xeon vCPUs at 2.10 GHz and 7 GB of RAM (comparable to a standard Colab runtime). Absolute timings depend on the hardware; the ratio between the algorithms is the stable quantity.
* Monte-Carlo errors are reported in the output (`MC_std_error`, replication counts in `Config`); use `--reps 2000` or more for publication-quality tables.

## 7. Troubleshooting

* *`ModuleNotFoundError`* - run `pip install -r requirements.txt` (and `matplotlib` for the base-paper figures).
* *`externally-managed-environment`* - add `--break-system-packages` or use a virtual environment.
* *Slow convergence experiment* - it loops over `k` up to 1000 and over the replications; use `--quick` or reduce `--reps`.
* *Excel chart not shown* - the chart is added with openpyxl; open the file in Excel or LibreOffice (not a text viewer).

## 8. Citation

If you use this code, please cite the article (preprint) and the base paper:

```
Al-Mohamad, M. Y., Parchami, A., Hydari, A. A. Quantile representation of LR random fuzzy numbers. Manuscript.
Parchami, A., Grzegorzewski, P., Romaniuk, M. (2024). Statistical simulations with LR random fuzzy numbers. Statistical Papers 65, 3585-3609.
```

## 9. Contact

Mohamad Yamen Al-Mohamad, Department of Statistics, University of Tabriz, Tabriz, Iran - yamenmohamad@tabrizu.ac.ir.
Corresponding author of the article: Abbas Parchami (Shahid Bahonar University of Kerman).
