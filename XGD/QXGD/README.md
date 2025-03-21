# The Quasi XGamma Distribution 
This document explains the implementation  of the Quasi XGamma Distribution (QXGD).
---

## Modules and Libraries Used

1. **NumPy**: For numerical operations like generating random numbers and performing mathematical computations.
2. **SciPy**: For statistical functions and optimization routines.
   - `scipy.stats` is used for distribution-related operations.
   - `scipy.optimize.minimize` is employed to maximize the log-likelihood function.
3. **Matplotlib**: For plotting the distributions of estimated parameters.

---

## Functions

### `dQXGD(x, alpha, theta)`
- **Purpose**: Computes the Probability Density Function (PDF) of the Quasi XGamma Distribution.
- **Parameters**:
  - `x`: A vector of values where the PDF is evaluated.
  - `alpha`, `theta`: Parameters of the Quasi XGamma distribution.
- **Output**: A vector representing the PDF values for each element of `x`.

### `pQXGD(x, alpha, theta)`
- **Purpose**: Computes the Cumulative Distribution Function (CDF) of the Quasi XGamma Distribution.
- **Parameters**:
  - `x`: A vector of values where the CDF is evaluated.
  - `alpha`, `theta`: Parameters of the Quasi XGamma distribution.
- **Output**: A vector representing the CDF values for each element of `x`.

### `rQXGD(n, alpha, theta)`
- **Purpose**: Generates a random sample from the Quasi XGamma Distribution using the inverse transform sampling method.
- **Parameters**:
  - `n`: Number of random samples to generate.
  - `alpha`, `theta`: Parameters of the Quasi XGamma distribution.
- **Output**: A vector of `n` random samples from the distribution.

### `log_likelihood(params, x)`
- **Purpose**: Computes the log-likelihood for the given parameters (α, θ) based on a sample `x` using the PDF of the Quasi XGamma Distribution.
- **Parameters**:
  - `params`: Parameters [α, θ] to be estimated.
  - `x`: A sample of data to compute the likelihood.
- **Output**: The negative log-likelihood value.

### `maximize_log_likelihood(x)`
- **Purpose**: Uses the `scipy.optimize.minimize` function to maximize the log-likelihood and estimate the parameters α and θ.
- **Parameters**:
  - `x`: A sample of data for which parameter estimation is performed.
- **Output**: The estimated parameters [α, θ] that maximize the likelihood.

### `monte_carlo_simulation(n, alpha_true, theta_true, num_simulations=1000)`
- **Purpose**: Performs the Monte Carlo simulation by generating `num_simulations` random samples and estimating the parameters using MLE. It then computes the bias and mean square error (MSE) of the estimates.
- **Parameters**:
  - `n`: Sample size.
  - `alpha_true`, `theta_true`: True values of the parameters for comparison.
  - `num_simulations`: Number of Monte Carlo simulations.
- **Output**: A dictionary containing:
    - `alpha_estimates`, `theta_estimates`: Estimated values for α and θ.
    - `alpha_bias`, `theta_bias`: Bias of the parameter estimates.
    - `alpha_mse`, `theta_mse`: Mean square error of the parameter estimates.

---

## Monte Carlo Simulation Setup

1. **Sample Size (`n`)**: The sample size for each simulation is set to 20.
2. **True Parameters**:
   - `alpha_true = 0.5`: The true value for the α parameter.
   - `theta_true = 0.5`: The true value for the θ parameter.
3. **Number of Simulations**: The number of Monte Carlo simulations is set to 10,000.

---

## Results and Interpretation

- **Bias and MSE**:
  - The bias and MSE values for the estimates of both α and θ are computed by comparing the true values with the estimated ones across all simulations.
  - A low bias and MSE indicate that the MLE method provides accurate and efficient parameter estimates.

```python
# Example Output:
Alpha Bias: 0.0001, Alpha MSE: 0.0102
Theta Bias: -0.0003, Theta MSE: 0.0115
