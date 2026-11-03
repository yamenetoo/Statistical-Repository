# Else

Two stand-alone notebooks that do not belong to the main experiment folders of this repository. They are independent of each other (different language, different topic). Each one runs top to bottom, is seeded for reproducibility, and keeps its saved outputs, so it can be read on GitHub without running anything.

| Notebook | Language | Topic |
|---|---|---|
| [`DLE1.ipynb`](DLE1.ipynb) | Python (TensorFlow/Keras) | Statistical pruning of a one-hidden-layer neural network: remove dead units, remove correlated units, transfer the weights, fine-tune |
| [`Dr homai.ipynb`](Dr%20homai.ipynb) | R | Density of `Z = W·X1 + (1−W)·X2` (normals mixed with a uniform weight): polynomial moment matching, maximum-likelihood normal, and the exact result |

## Contents

1. [DLE1.ipynb](#1-dle1ipynb-statistical-pruning-of-a-neural-network)
2. [Dr homai.ipynb](#2-dr-homaiipynb-density-of-a-normaluniform-mixture)
3. [How to run](#3-how-to-run)
4. [Requirements](#4-requirements)
5. [Changes from the first version](#5-changes-from-the-first-version)

---

## 1. `DLE1.ipynb`: statistical pruning of a neural network

**Idea.** A one-hidden-layer binary classifier is a logistic regression whose regressors are the hidden activations. Treat those activations as a data matrix (samples × units) and use descriptive statistics to decide which units are unnecessary.

**Data and baseline.** Iris restricted to setosa and versicolor (100 samples, 4 features), split 80 / 20 with `random_state=42`. Baseline model: `Input(4) → Dense(16, ReLU, "HL") → Dense(1, sigmoid, "OL")`, 97 trainable parameters, Adam, binary cross-entropy, 10 epochs, batch size 1.

### Structure

| Section | Content |
|---|---|
| 1. Setup | Imports, configuration constants (seed, epochs, thresholds), random seeds |
| 2. Helper functions | `build_model`, `hidden_activations`, `find_dead_units`, `select_uncorrelated_units`, `prune_hidden_units`, `summarise`, plotting helpers, optional architecture diagram |
| 3. Data and baseline model | Load data, split, train the 16-unit network |
| 4. Hidden-layer diagnostics | Box plots, descriptive statistics and correlation matrix of the 100 × 16 activation matrix |
| 5. Remove dead units | Units with zero activation variance are deleted together with their bias and outgoing weight; predictions are checked against the baseline |
| 6. Remove correlated units | Greedy filter (a unit is kept only if its &#124;r&#124; with every kept unit is at most 0.7), weight transfer, scatter plot and single-unit logistic regression |
| 7. Fine-tune | 25 more epochs on the pruned network |
| 8. Comparison | One table with units, parameters, test loss and test accuracy for every stage |

`prune_hidden_units(model, keep)` copies, for each kept unit, its column of the first kernel, its bias, and its row of the output kernel, using 0-based indices throughout.

### Results (seed 42, as stored in the notebook)

| Stage | Hidden units | Parameters | Test loss | Test accuracy |
|---|---:|---:|---:|---:|
| Baseline | 16 | 97 | 0.0322 | 1.00 |
| Dead units removed | 6 | 37 | 0.0322 | 1.00 |
| Correlated units removed, before fine-tuning | 1 | 7 | 0.5721 | 0.50 |
| Pruned and fine-tuned | 1 | 7 | 0.2244 | 1.00 |

- 10 of the 16 hidden units were dead (units 0, 2, 3, 4, 7, 9, 10, 11, 13, 14). Removing them is lossless: the largest prediction difference on the test set is 1.5e-8.
- The six remaining units are strongly correlated ((&#124;r&#124; up to 0.99)#124;r(&#124;r&#124; up to 0.99)#124; close to 1 for most pairs), so the filter keeps a single unit. Without retraining the function changes (accuracy 0.50), and 25 epochs of fine-tuning recover test accuracy 1.00 with 7 parameters instead of 97.
- That single unit separates the two classes perfectly, so the maximum-likelihood logistic regression on it does not exist (the notebook reports this instead of failing).

The exact numbers depend on the seed and on the TensorFlow version; the qualitative picture (many dead units, strong redundancy among the rest) is what the notebook demonstrates.

---

## 2. `Dr homai.ipynb`: density of a normal–uniform mixture

**Question.** Let `X1, X2 ~ N(0,1)` and `W ~ U(0,1)` be independent. What is the density of `Z = W·X1 + (1−W)·X2`? Conditional on `W = w`, `Z ~ N(0, w² + (1−w)²)`, so `Z` is a scale mixture of normals.

The notebook simulates `n = 1000` observations (`set.seed(123)`) and compares:

| Section | Method | Details |
|---|---|---|
| 3 | Polynomial moment matching | Degree-10 polynomial on [−3, 3] whose moments equal the sample moments of `Z`; solves `M w = m` with `M[r,j] = ∫ x^(r+j) dx`. Reports the condition number, the minimum of the fit and its mass |
| 4 | Parametric fit | Normal distribution by maximum likelihood (`fitdistrplus::fitdist`) |
| 5 | Theory | `E(Z) = 0`, `Var(Z) = 2/3`, `E(Z⁴) = 7/5` (kurtosis 3.15, excess 0.15), and the exact density `f_Z(z) = ∫₀¹ φ(z; sd = √(w² + (1−w)²)) dw` by numerical integration |
| 6 | Comparison | Plot of all curves over the histogram, and the integrated squared error (ISE) of each against the exact density |

### Results (as stored in the notebook)

| Quantity | Theory | Sample |
|---|---:|---:|
| Mean | 0 | 0.0249 |
| Variance | 0.6667 | 0.7151 |
| Excess kurtosis | 0.15 | 0.0896 |

| Approximation | ISE vs exact density |
|---|---:|
| N(0, 2/3), same mean and variance | 7.8e-05 |
| Normal, maximum likelihood (sd 0.845) | 8.0e-04 |
| Polynomial, N = 10 | 5.4e-04 |

The polynomial system has condition number 2.3e10, and the fitted polynomial dips to −0.019 (it is not a valid density), with mass 0.9999 on [−3, 3].

---

## 3. How to run

**Google Colab.** Open the notebook on GitHub and use the Colab badge in its first cell. Choose the runtime that matches the notebook: Python 3 for `DLE1`, R for `Dr homai` (Runtime → Change runtime type).

**Locally.**

```bash
git clone https://github.com/yamenetoo/Statistical-Repository.git
cd Statistical-Repository/Else
jupyter lab
```

The R notebook needs the R kernel: `install.packages("IRkernel"); IRkernel::installspec()`.

## 4. Requirements

| Notebook | Packages |
|---|---|
| `DLE1.ipynb` | Python 3, `tensorflow`/`keras`, `numpy`, `pandas`, `matplotlib`, `scikit-learn`, `statsmodels`; optional `plot-ann` for the architecture diagrams (`pip install plot-ann`) |
| `Dr homai.ipynb` | R, `ggplot2`, `fitdistrplus` (the install line is commented out in the first cell) |

## 5. Changes from the first version

Both notebooks were rewritten as structured, seeded, function-based analyses and re-executed, so all stored outputs come from the current code. Substantive corrections:

**`DLE1.ipynb`**

- Pruning is now done by one function that works on explicit 0-based unit indices. The first version mixed 1-based column labels with 0-based indices when building the third model, so it kept the wrong units.
- Model accuracies are computed for every stage with the same helper; the first version printed Model 2's accuracy in the Model 3 cell and never evaluated Model 3.
- The correlation filter is a documented greedy rule (keep a unit only if it is not strongly correlated with an already kept one), and dead units use a numerical tolerance on the standard deviation instead of reading a fixed row of `describe()`.
- The unused variables, the global-variable dependency in the box-plot helper, and the plotting inside the model builder were removed; seeds were added.

**`Dr homai.ipynb`**

- The first version concluded that `Z` is exactly `N(0, 2/3)`. The mean and variance are right, but `Z` is a scale mixture with excess kurtosis 0.15, so the conclusion was changed and the exact density is now computed.
- The "exact" normal curve called `dnorm(sd = 2/3)`, but 2/3 is the variance; the standard deviation is `√(2/3) ≈ 0.816`. This is corrected.
- The polynomial fit now reports its condition number, its minimum and its mass; the unused spline step and the empty cell were removed.

The Colab badges now point to this repository.
