# Else

Two stand-alone teaching/research notebooks that do not belong to the main experiment folders of this repository. They are independent of each other (different language, different topic), and each one runs top to bottom in Google Colab.

| Notebook | Language | Topic | Size |
|---|---|---|---|
| [`DLE1.ipynb`](DLE1.ipynb) | Python (TensorFlow/Keras) | Statistical pruning of a one-hidden-layer neural network: remove dead neurons, remove highly correlated neurons, transfer the weights | 62 cells, about 700 KB |
| [`Dr homai.ipynb`](Dr%20homai.ipynb) | R | Approximating the density of `Z = W·X1 + (1−W)·X2` (normals mixed with a uniform weight) by a moment-matching polynomial, by `fitdistrplus`, and by theory | 18 cells, about 290 KB |

Both notebooks keep their saved outputs (tables, plots, training logs), so they can be read on GitHub without running anything.

---

## Contents

1. [DLE1.ipynb: shrinking a network with statistics](#1-dle1ipynb-shrinking-a-network-with-statistics)
2. [Dr homai.ipynb: density of a normal–uniform mixture](#2-dr-homaiipynb-density-of-a-normaluniform-mixture)
3. [How to run](#3-how-to-run)
4. [Requirements](#4-requirements)
5. [Notes and known issues](#5-notes-and-known-issues)

---

## 1. `DLE1.ipynb`: shrinking a network with statistics

**Idea.** Treat the activations of the hidden layer as a data matrix (samples × neurons) and use ordinary descriptive statistics on it to decide which neurons are unnecessary. The output layer of the small binary classifier is a logistic regression whose independent variables are those activations, so removing a neuron is the same as removing a regressor.

### Data and base model

- Data: Iris, restricted to classes 0 and 1 (setosa and versicolor), 100 samples, 4 features.
- Split: 80 train / 20 test, `random_state=42`.
- Model 1: `Input(4) → Dense(16, relu, "HL") → Dense(1, sigmoid, "OL")`, 97 parameters, Adam, binary cross-entropy, 10 epochs, batch size 1.
- Result: test accuracy 1.0; training loss goes from 0.76 to 0.18 over the 10 epochs.

### Pipeline

| Step | What the notebook does | Result in the saved outputs |
|---|---|---|
| 1. Helper functions | `drop_node` (delete units from a weight array), `intermediate_layer` (activations of a named layer), `plot_boxplot`, `plot_corr_matrix`, `plot_scatter`, `generate_pairs`, `plot_ann` (network diagram via the `plot-ann` package), `make_model` | |
| 2. Inspect the hidden layer | Activations of `HL` on all 100 samples: box plots, `describe()`, correlation heat map | 100 × 16 activation matrix |
| 3. Find dead neurons | Neurons whose standard deviation is exactly 0 on the data (always output 0 after ReLU) | 6 of 16: units 2, 3, 8, 9, 11, 14 |
| 4. Build Model 2 | New network with 16 − 6 = 10 hidden units. Weights are copied from Model 1 after deleting the dead units' input columns, their biases and their output rows | Model 2 predictions equal Model 1 to about 1e-7; test accuracy 1.0. The removal is lossless |
| 5. Correlated neurons | Pairwise correlation of the 10 remaining activations; for each pair with &#124;r&#124; > 0.7 the later column is dropped | 3 units remain (columns 1, 2 and 7 of Model 2, counted from 1) |
| 6. Look at the survivors | Scatter plots of the remaining pairs (the classes are completely separated), and a logistic regression (`statsmodels.Logit`) on the unit with the largest variance | Coefficient −1.3543 (std. err. 0.301, z = −4.50, p < 0.001), pseudo R² = 0.3804, no constant term |
| 7. Build Model 3 | 3 hidden units, weights transferred from Model 2 | Immediately after the transfer the predictions change (for example 0.96 becomes 0.29 for the first test row) |
| 8. Fine-tune Model 3 | 25 more epochs, batch size 1 | Training accuracy 1.0, loss 0.033. The printed test predictions are close to 1 (at least 0.96) for class-1 rows and close to 0 (at most 0.06) for class-0 rows |

Parameter counts (from the layer sizes): Model 1 has 97, Model 2 has 61, Model 3 has 19.

### What the notebook shows

- Dead ReLU units can be removed exactly, with no retraining.
- Strongly correlated hidden units carry largely redundant information; dropping them changes the function, but a short fine-tune recovers it, here with a network about five times smaller than the original.

---

## 2. `Dr homai.ipynb`: density of a normal–uniform mixture

**Question.** Let `X1, X2 ~ N(0,1)` independent and `W ~ U(0,1)` independent of both. What is the density of

```
Z = W·X1 + (1 − W)·X2 ?
```

The notebook (R kernel) generates `n = 1000` observations with `set.seed(123)` and compares three descriptions of the density of `Z`.

| Section | Method | Details |
|---|---|---|
| Data | Simulation | `X1`, `X2`, `W`, `Z` as above; histograms with kernel density estimates for `X1`, `X2`, `Z` |
| 1 | Polynomial moment matching | Density approximated on [−3, 3] by a degree N = 10 polynomial `f(x) = Σ w_j x^j` whose raw moments match the sample moments of `Z`. The weights solve `M w = m` with `M[i,j] = (b^(i+j−1) − a^(i+j−1)) / (i+j−1)` (integrals of monomials) and `m_r` = sample mean of `Z^r`, r = 0..10. The curve is smoothed with `spline` and drawn over the histogram |
| 2 | Parametric fit | `fitdistrplus::fitdist(Z, "norm")`: maximum-likelihood normal fit; its density is drawn next to the polynomial fit |
| 3 | Theory | `E(Z) = 0` and `Var(Z) = E(W²)+E((1−W)²) = 2/3`; the theoretical normal curve is added to the comparison plot |

The notebook has an explanation in Markdown for each code block (what the code does line by line).

---

## 3. How to run

**Google Colab (simplest).** Open the notebook on GitHub and use the Colab badge in its first cell. Choose the runtime that matches the notebook: Python 3 for `DLE1`, R for `Dr homai` (Runtime → Change runtime type).

**Locally.**

```bash
git clone https://github.com/yamenetoo/Statistical-Repository.git
cd Statistical-Repository/Else
jupyter lab            # or: jupyter notebook
```

For the R notebook install the R kernel once: `install.packages("IRkernel"); IRkernel::installspec()`.

---

## 4. Requirements

| Notebook | Packages |
|---|---|
| `DLE1.ipynb` | Python 3.10 (as stored in the notebook), `tensorflow`/`keras`, `numpy`, `pandas`, `scikit-learn`, `matplotlib`, `statsmodels`, `plot-ann` (installed in the first cell with `!pip install plot-ann`) |
| `Dr homai.ipynb` | R, `fitdistrplus`, `ggplot2` (the install line is commented out in the notebook; uncomment it on the first run) |

The saved outputs were produced in Colab in November 2023. A newer TensorFlow/Keras version may print slightly different training logs.

---

## 5. Notes and known issues

These were found while reading the notebooks. The notebooks themselves were not changed.

**`Dr homai.ipynb`**

- **`Z` is not exactly normal.** Conditional on `W`, `Z ~ N(0, W² + (1−W)²)`, so `Z` is a scale mixture of normals. The mean (0) and the variance (2/3) in the theory section are correct, but the conclusion that `Z` is "a normal distribution with mean 0 and variance 2/3" is not: `E(Z⁴) = 3·E[(W²+(1−W)²)²] = 7/5`, so the kurtosis is 3.15 (excess kurtosis 0.15), not 3. A simulation of 2 million draws gives variance 0.667 and excess kurtosis 0.150. The normal curve is a close but not an exact description.
- **Wrong scale in the "exact" curve.** The theory section calls `dnorm(x_values, mean = 0, sd = 2/3)`, but 2/3 is the variance. The standard deviation is `sqrt(2/3) ≈ 0.816`. The plotted curve labelled `N(0,2/3)` is therefore narrower than it should be. The `fitdist` curve (sd estimated from the data, about 0.82) is the correct one.
- The polynomial method solves a Hilbert-type linear system, which is ill-conditioned for N = 10, and the fitted polynomial is not guaranteed to be non-negative.
- The last code cell is empty.

**`DLE1.ipynb`**

- **Off-by-one when building Model 3.** The kept columns are labelled 1, 2, 7 (counted from 1), but the code that builds `empty_col` compares those labels with indices counted from 0. As a result, Model 3 keeps the units at 0-based positions 1, 2 and 7 instead of 0, 1 and 6. This probably explains why the predictions change so much right after the transfer, before the fine-tuning.
- **Wrong variable in a print.** After building Model 3, the accuracy line prints `accuracy2` (Model 2), so the `1.0` shown there is not Model 3's accuracy. Model 3's test accuracy after fine-tuning is never computed; only its predicted probabilities are shown.
- **Unused or misplaced code.** `dic3` is created but the loop writes into `dic`; the column `Out` built with `sm.add_constant` is not used by the logistic regression, which is fitted without an intercept; `plot_boxplot` takes its number of labels from the global `layer_output_2`.
- **A near-dead unit survives.** Column 7 of Model 2 is zero for 99 of the 100 samples (its largest value is 0.0228) but was not removed as a dead unit because its standard deviation is not exactly 0; it passes the correlation filter because its correlation with the other units is small.
- Dead-neuron detection reads `summary_stats.iloc[2]`, which is the standard-deviation row of `describe()`; it silently depends on that row order.
- The Colab badges at the top of both notebooks point to a different repository (`yamenetoo/ST`), not to this one, so they open copies stored there.
- Some labels and comments are in Arabic (a scatter-plot legend and one cell comment); a plot font with Arabic support is needed to render them.
