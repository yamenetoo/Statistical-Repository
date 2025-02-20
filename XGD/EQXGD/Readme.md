# Exponentiated Quasi X-Gamma Distribution (EQXGD)

This Python module implements functions to work with the **Exponentiated Quasi X-Gamma Distribution (EQXGD)**. It provides the following functionality:

- **Cumulative Distribution Function (CDF)**: Evaluates the CDF for given input values.
- **Probability Density Function (PDF)**: Evaluates the PDF for given input values.
- **Random Variable Generation**: Generates random variables following the EQXGD distribution.
- **Maximum Likelihood Estimation (MLE)**: Fits parameters of the EQXGD distribution to data using MLE.
- **Monte Carlo Simulation**: Simulates MLE fitting for multiple simulations to estimate the mean and standard deviation of the parameters.
- **Fit and Plot**: Fits the EQXGD distribution to a data sample and visualizes the result using a histogram and the fitted density.

## Functions

### `pEQXGD(x, theta, alpha, beta)`
Evaluates the **Cumulative Distribution Function (CDF)** for the **Exponentiated Quasi X-Gamma Distribution (EQXGD)**.

#### Parameters:
- `x` : array-like  
  Values at which to evaluate the CDF.
- `theta` : float  
  Scale parameter (θ > 0).
- `alpha` : float  
  Shape parameter (α > 0).
- `beta` : float  
  Exponentiation parameter (β > 0).

#### Returns:
- `array-like`  
  CDF values at the specified points `x`.

### `dEQXGD(x, theta, alpha, beta)`
Evaluates the **Probability Density Function (PDF)** for the **Exponentiated Quasi X-Gamma Distribution (EQXGD)**.

#### Parameters:
- `x` : array-like  
  Values at which to evaluate the PDF.
- `theta` : float  
  Scale parameter (θ > 0).
- `alpha` : float  
  Shape parameter (α > 0).
- `beta` : float  
  Exponentiation parameter (β > 0).

#### Returns:
- `array-like`  
  PDF values at the specified points `x`.

### `rEQXGD(theta, alpha, beta, n)`
Generates `n` random variables following the **Exponentiated Quasi X-Gamma Distribution (EQXGD)**.

#### Parameters:
- `theta` : float  
  Scale parameter (θ > 0).
- `alpha` : float  
  Shape parameter (α > 0).
- `beta` : float  
  Exponentiation parameter (β > 0).
- `n` : int  
  The number of random variables to generate.

#### Returns:
- `array-like`  
  `n` random variables following the EQXGD distribution.

### `neg_log_likelihood(params, samples)`
Computes the **Negative Log-Likelihood** for the **Exponentiated Quasi X-Gamma Distribution (EQXGD)**, used for parameter estimation.

#### Parameters:
- `params` : list of floats  
  A list of parameters `[theta, alpha, beta]` to optimize.
- `samples` : array-like  
  The data to fit the distribution to.

#### Returns:
- float  
  The negative log-likelihood value.

### `monte_carlo_simulation(true_params, sample_size, num_simulations)`
Performs a **Monte Carlo simulation** to estimate the **Maximum Likelihood Estimation (MLE)** of the parameters of the EQXGD distribution across multiple simulations.

#### Parameters:
- `true_params` : tuple  
  The true parameters `[theta_true, alpha_true, beta_true]` used to generate samples.
- `sample_size` : int  
  The number of samples to generate in each simulation.
- `num_simulations` : int  
  The number of simulations to run.

#### Returns:
- dict  
  Estimated mean and standard deviation of each parameter (`theta`, `alpha`, `beta`).

### `fitAndPlot(data, bins=None)`
Fits the **Exponentiated Quasi X-Gamma Distribution (EQXGD)** to a data sample using MLE, then visualizes the data and fitted density.

#### Parameters:
- `data` : array-like  
  The data sample to fit.
- `bins` : int, optional  
  The number of bins to use for the histogram.

#### Returns:
- None  
  Displays the plot and prints the MLE parameters, log-likelihood, AIC, and BIC values.

## Dependencies

This module requires the following Python packages:
- `numpy`
- `scipy`
- `matplotlib` (for plotting in `fitAndPlot` function)

You can install the required dependencies using `pip`:

```bash
pip install numpy scipy matplotlib
```
  
