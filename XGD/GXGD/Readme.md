 


# GXGD Distribution

The **GXGD (Generalized XGamma- Distribution)** is a flexible probability distribution that arises as a mixture of an **exponential distribution** and a **gamma distribution**, weighted by parameters derived from a geometric mixing mechanism.

## Probability Density Function (PDF)

The PDF of the GXGD distribution is defined as:

```
f_X(x; θ, a, b, w) = (w / (1 + w)) * θ * e^(-θ * x) 
                   + (1 / (1 + w)) * (b^a / Γ(a)) * x^(a-1) * e^(-b * x),
                   where x > 0, θ > 0, a > 0, b > 0, w > 0.
```

Where:
- `(w / (1 + w))` and `(1 / (1 + w))` are the mixing weights for the exponential and gamma components, respectively.
- `θ * e^(-θ * x)` is the PDF of the exponential distribution.
- `(b^a / Γ(a)) * x^(a-1) * e^(-b * x)` is the PDF of the gamma distribution.
- `Γ(a)` is the complete gamma function.

## Cumulative Distribution Function (CDF)

The CDF of the GXGD distribution is given by:

```
F_X(x; θ, a, b, w) = (w / (1 + w)) * (1 - e^(-θ * x)) 
                   + (1 / (1 + w)) * (γ(a, b * x) / Γ(a)),
                   where x > 0.
```

Where:
- `γ(a, b * x)` is the lower incomplete gamma function.
- `Γ(a)` is the complete gamma function.

## Key Properties

1. **Support**:  
   The GXGD distribution is defined for `x > 0`, making it suitable for modeling positive-valued data such as survival times, waiting times, or other non-negative random variables.

2. **Mixing Mechanism**:  
   The parameter `w > 0` controls the relative contribution of the exponential and gamma components:
   - When `w → 0`, the GXGD distribution approaches the gamma distribution (`X ~ Gamma(a, b)`).
   - When `w → ∞`, the GXGD distribution approaches the exponential distribution (`X ~ Exp(θ)`).

3. **Moments**:  
   The moments of the GXGD distribution can be derived as weighted sums of the moments of the exponential and gamma components:
   - Mean:  
     ```
     E[X] = (w / (1 + w)) * (1 / θ) + (1 / (1 + w)) * (a / b)
     ```
   - Variance:  
     ```
     Var(X) = (w / (1 + w)) * (1 / θ^2) + (1 / (1 + w)) * (a / b^2)
     ```

4. **Flexibility**:  
   The GXGD distribution combines the simplicity of the exponential distribution with the flexibility of the gamma distribution, allowing it to model a wide range of data shapes, including right-skewed distributions commonly observed in survival analysis and reliability studies.

## Applications

The GXGD distribution is particularly useful in scenarios where the data exhibit characteristics of both exponential and gamma distributions, such as:
- **Survival Analysis**: Modeling time-to-event data where early failures may follow an exponential pattern, while later failures follow a gamma pattern.
- **Queueing Theory**: Modeling service times or inter-arrival times with varying degrees of variability.
- **Reliability Engineering**: Capturing systems with mixed failure modes.

## Parameters

- `θ > 0`: The rate parameter of the exponential component.
- `a > 0`: The shape parameter of the gamma component.
- `b > 0`: The rate parameter of the gamma component.
- `w > 0`: The mixing weight parameter controlling the balance between the exponential and gamma components.

## Usage

The GXGD distribution can be used for:
- Parameter estimation using Maximum Likelihood Estimation (MLE).
  ```
  data=[]
  initial_params=[0.5, 1.0, 1.5, 1.0]
  learning_rate=0.01,
  epochs=1000,
  fit_params = SGD(data, initial_params, learning_rate=learning_rate, epochs=epochs, ver=False)
  ```
- Monte Carlo simulations to evaluate the performance of estimators.
 ```
theta_true = 0.5
a_true = 1.5
b_true = 1.5
w_true = 2.5
results = monte_carlo_simulation(theta_true, a_true, b_true, w_true,n_simulations=100,
                                  n_samples_per_simulation=50,initial_params=[0.5, 1.0, 1.5, 1.0],
                                  learning_rate=0.01, epochs=1000,
                                  verbose=False)
print("\nMonte Carlo Simulation Results:")
print(f"True Parameters: θ={theta_true}, a={a_true}, b={b_true}, w={w_true}")
print(f"Estimated θ: Mean={results['theta']['mean']:.4f}, Std={results['theta']['std']:.4f}")
print(f"Estimated a: Mean={results['a']['mean']:.4f}, Std={results['a']['std']:.4f}")
print(f"Estimated b: Mean={results['b']['mean']:.4f}, Std={results['b']['std']:.4f}")
print(f"Estimated w: Mean={results['w']['mean']:.4f}, Std={results['w']['std']:.4f}")
 ```
- Goodness-of-fit testing for real-world datasets.

## References
 
