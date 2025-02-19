import numpy as np
import scipy.stats as stats
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Probability Density Function (PDF) of Quasi XGamma Distribution
def dQXGD(x, alpha, theta):
    if np.any(x <= 0) or alpha <= 0 or theta <= 0:
        return np.zeros_like(x)
    
    numerator = theta * (alpha + (theta**2 * x**2) / 2)
    denominator = 1 + alpha
    exp_term = np.exp(-theta * x)
    
    return (numerator / denominator) * exp_term

# Cumulative Distribution Function (CDF) of Quasi XGamma Distribution
def pQXGD(x, alpha, theta):
    if np.any(x <= 0) or alpha <= 0 or theta <= 0:
        return np.zeros_like(x)
    
    numerator = 1 + alpha + theta * x + (theta**2 * x**2) / 2
    denominator = 1 + alpha
    exp_term = np.exp(-theta * x)
    
    return 1 - (numerator / denominator) * exp_term

# Random Sample Generator from QXGD
def rQXGD(n, alpha, theta):
    U = np.random.uniform(0, 1, n)  # Step 1
    V = np.random.exponential(1/theta, n)  # Step 2
    W = np.random.gamma(3, 1/theta, n)  # Step 3
    X = np.where(U <= alpha / (1 + alpha), V, W)  # Step 4
    return X

# Log-Likelihood Function
def log_likelihood(params, x):
    alpha, theta = params
    if alpha <= 0 or theta <= 0:
        return np.inf  # Ensure positive parameters
    return -np.sum(np.log(dQXGD(x, alpha, theta)))  # Negative Log-Likelihood

# MLE Estimation using Optimization
def maximize_log_likelihood(x):
    result = minimize(log_likelihood, [0.5, 0.5], args=(x,), method="Nelder-Mead")
    return result.x if result.success else None

# Monte Carlo Simulation
def monte_carlo_simulation(n, alpha_true, theta_true, num_simulations=1000):
    alpha_estimates = []
    theta_estimates = []
    
    for _ in range(num_simulations):
        sample = rQXGD(n, alpha_true, theta_true)  # Generate sample
        estimates = maximize_log_likelihood(sample)  # Estimate parameters
        
        if estimates is not None:
            alpha_estimates.append(estimates[0])
            theta_estimates.append(estimates[1])

    # Compute Bias and MSE
    alpha_bias = np.mean(alpha_estimates) - alpha_true
    theta_bias = np.mean(theta_estimates) - theta_true
    alpha_mse = np.mean((np.array(alpha_estimates) - alpha_true) ** 2)
    theta_mse = np.mean((np.array(theta_estimates) - theta_true) ** 2)

    return {
        "alpha_estimates": alpha_estimates,
        "theta_estimates": theta_estimates,
        "alpha_bias": alpha_bias,
        "theta_bias": theta_bias,
        "alpha_mse": alpha_mse,
        "theta_mse": theta_mse
    }

# Simulation Parameters
n = 20  # Sample Size
alpha_true = .5  # True Alpha
theta_true = .5 # True Theta
num_simulations = 10000  # Number of Simulations

# Run Monte Carlo Simulation
results = monte_carlo_simulation(n, alpha_true, theta_true, num_simulations)

# Print Bias and MSE
print(f"Alpha Bias: {results['alpha_bias']:.4f}, Alpha MSE: {results['alpha_mse']:.4f}")
print(f"Theta Bias: {results['theta_bias']:.4f}, Theta MSE: {results['theta_mse']:.4f}")

# Plot Distribution of Estimates
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.hist(results["alpha_estimates"], bins=30, alpha=0.7, color='blue', edgecolor='black')
plt.axvline(alpha_true, color='red', linestyle='dashed', linewidth=2, label=f'True α={alpha_true}')
plt.xlabel("Alpha Estimates")
plt.ylabel("Frequency")
plt.title("Distribution of Alpha Estimates")
plt.legend()

plt.subplot(1, 2, 2)
plt.hist(results["theta_estimates"], bins=30, alpha=0.7, color='green', edgecolor='black')
plt.axvline(theta_true, color='red', linestyle='dashed', linewidth=2, label=f'True θ={theta_true}')
plt.xlabel("Theta Estimates")
plt.ylabel("Frequency")
plt.title("Distribution of Theta Estimates")
plt.legend()
plt.tight_layout()
plt.show()
