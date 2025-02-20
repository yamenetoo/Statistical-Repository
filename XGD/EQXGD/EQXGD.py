import numpy as np
from scipy.optimize import brentq
from scipy.optimize import minimize

def pEQXGD(x, theta, alpha, beta):
    """
    Cumulative Distribution Function (CDF) of Exponentiated Quasi X-Gamma Distribution.
    Parameters:
    x : array-like
        Values at which to evaluate the CDF.
    theta : float
        Scale parameter (theta > 0).
    alpha : float
        Shape parameter (alpha > 0).
    beta : float
        Exponentiation parameter (beta > 0).
    Returns:
    array-like : CDF values
    """
    if np.any(x <= 0) or theta <= 0 or alpha <= 0 or beta <= 0:
        return np.zeros_like(x)
    return np.power(pQXGD(x,alpha=alpha,theta=theta), beta)
def dEQXGD(x, theta, alpha, beta):
    """
    Probability Density Function (PDF) of Exponentiated Quasi X-Gamma Distribution.

    Parameters:
    x : array-like
        Values at which to evaluate the PDF.
    theta : float
        Scale parameter (theta > 0).
    alpha : float
        Shape parameter (alpha > 0).
    beta : float
        Exponentiation parameter (beta > 0).

    Returns:
    array-like : PDF values
    """
    if np.any(x <= 0) or theta <= 0 or alpha <= 0 or beta <= 0:
        return np.zeros_like(x)
    F_QXGD = pQXGD(x,alpha=alpha,theta=theta)
    f_QXGD = dQXGD(x,alpha=alpha,theta=theta)
    return beta * np.power(F_QXGD, beta - 1) * f_QXGD

def rEQXGD(theta, alpha, beta,n):
    """
    Quantile function (inverse CDF) of the Exponentiated Quasi X-Gamma Distribution (EQXGD).
    Solves for x numerically given p.
    """
    p=np.random.uniform(0,1,n)
    def root_func(x, p):
        return pEQXGD(x, theta, alpha, beta) - p
    lower_bound = 0.001
    upper_bound = 10
    result = []
    for pi in p:
        try:
            root = brentq(root_func, lower_bound, upper_bound, args=(pi))
            result.append(root)
        except ValueError as e:
            result.append(np.nan)  # Append NaN if there's an issue with root-finding
    return np.array(result)

def neg_log_likelihood(params, samples):
    theta, alpha, beta = params
    if theta <= 0 or alpha <= 0 or beta <= 0:
        return np.inf   
    n = len(samples)
    log_likelihood = 0
    for x in samples:
        F_QXGD = pQXGD(x, alpha=alpha, theta=theta)
        f_QXGD = dQXGD(x, alpha=alpha, theta=theta)
        if F_QXGD == 0 or f_QXGD == 0:
            return np.inf   
        log_likelihood += np.log(beta) + (beta - 1) * np.log(F_QXGD) + np.log(f_QXGD)
    return -log_likelihood

def monte_carlo_simulation(true_params, sample_size, num_simulations):
    theta_true, alpha_true, beta_true = true_params
    theta_mles = []
    alpha_mles = []
    beta_mles = []
    for _ in range(num_simulations):
        samples = rEQXGD(theta_true, alpha_true, beta_true, sample_size)
        initial_params = [0.5, 0.5, 0.5]
        result = minimize(neg_log_likelihood, initial_params, args=(samples,), 
                          method='L-BFGS-B', bounds=[(1e-6, None), (1e-6, None), (1e-6, None)])
        theta_mle, alpha_mle, beta_mle = result.x
        theta_mles.append(theta_mle)
        alpha_mles.append(alpha_mle)
        beta_mles.append(beta_mle)
    theta_mean = np.mean(theta_mles)
    alpha_mean = np.mean(alpha_mles)
    beta_mean = np.mean(beta_mles)
    theta_std = np.std(theta_mles)
    alpha_std = np.std(alpha_mles)
    beta_std = np.std(beta_mles)
    
    return {
        'theta': {'mean': theta_mean, 'std': theta_std},
        'alpha': {'mean': alpha_mean, 'std': alpha_std},
        'beta': {'mean': beta_mean, 'std': beta_std}
    }

def fitAndPlot(data,bins=None):
    samples=data
    initial_params = [0.5, 0.5, 0.5]
    result = minimize(neg_log_likelihood, initial_params, args=(data,),method='L-BFGS-B', bounds=[(1e-6, None), (1e-6, None), (1e-6, None)])
    theta_mle, alpha_mle, beta_mle = result.x
    hist=plt.hist(samples, density=True,bins=bins, alpha=0.6, color='b', edgecolor='black')
    x_vals = np.linspace(min(samples), max(samples) + 0.2, 100)
    plt.plot(x_vals, dEQXGD(x_vals,theta=theta_mle,alpha=alpha_mle,beta=0.9), 'r-', label='Density')
    plt.ylim(0, hist[0].max()+0.001)
    plt.xlim(min(hist[1]), max(hist[1]) + 0.2)
    plt.show()
    log_likelihood = -result.fun  # since we minimized the negative log-likelihood
    k = 3  # number of parameters
    n = len(samples)  # number of data points
    # AIC and BIC calculations
    AIC = 2 * k - 2 * log_likelihood
    BIC = np.log(n) * k - 2 * log_likelihood
    print(f"MLE Parameters: θ = {theta_mle}, α = {alpha_mle}, β = {beta_mle}")
    print(f"log_likelihood: {-log_likelihood}")
    print(f"AIC: {AIC}")
    print(f"BIC: {BIC}")
