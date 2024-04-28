# Bayesian Estimation for Topp-Leone Distribution

This repository contains R functions for computing Bayesian estimators under different loss functions and prior distributions for n for Topp-Leone Distribution. The implemented loss functions include:

1. **Pointwise Logarithmic Loss Function (PLF)**
2. **Weighted Loss Function (WLF)**
3. **Entropy Loss Function (ELF)**
4. **Squared-Log Error Loss Function (SLELF)**

These functions are designed to estimate parameters in a Bayesian setting under various prior distributions:

- **Uniform prior**
- **Jeffreys prior**
- **Exponential prior**
- **Gamma prior**
## NOTE
	Let $x_{(1)}, x_{(2)}, x_{(3)}, \ldots, x_{(r)}$ be the ordered observations. Then $$\xi_{ir} = -\left\{\sum_{i=1}^{r} \ln(2x_{(i)} - x_{(i)}^2) + k \ln(2x_{(r)} - x_{(r)}^2)\right\}$$.
	
	\textbf{PLF loss function}
	
	The Bayes estimator under uniform prior is:
	$
	\hat{V}_{plf\_U} = \sqrt{\frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+3)}{\xi_{ir}^{r+3}}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+1)}{\xi_{ir}^{r+1}}}}
	$
	
	The Bayes estimator under Jeffreys prior is:
	$
	\hat{V}_{\text{plf\_Jeffreys}} = \sqrt{\frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+2)}{\xi_{ir}^{r+2}}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r)}{\xi_{ir}^r}}}
	$
	
	The Bayes estimator under exponential prior is:
	$
	\hat{V}_{\text{PLF\_EXP}} = \sqrt{\frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+3)}{(w + \xi_{ir})^{r+3}}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+1)}{(w + \xi_{ir})^{r+1}}}}
	$
	
	The Bayes estimator under gamma prior is:
	$
	\hat{V}_{\text{PLF\_gamma}} = \sqrt{\frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+a+2)}{(b + \xi_{ir})^{r+a+2}}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+a)}{(b + \xi_{ir})^{r+a}}}}
	$
	
	\textbf{Weighted Loss Function (WLF)}
	
	The Bayes estimator under uniform prior is:
	$
	\hat{V}_{\text{WLF\_UNI}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+1)}{\xi_{ir}^{r+1}}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r)}{\xi_{ir}^r}}
	$
	
	The Bayes estimator under Jeffreys prior is:
	$
	\hat{V}_{\text{WLF\_Jeffreys}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r)}{\xi_{ir}^r}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r-1)}{\xi_{ir}^{r-1}}}
	$
	
	The Bayes estimator under exponential prior is:
	$
	\hat{V}_{\text{WLF\_EXP}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+1)}{(w + \xi_{ir})^{r+1}}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r)}{(w + \xi_{ir})^r}}
	$
	
	The Bayes estimator under gamma prior is:
	$
	\hat{V}_{\text{WLF\_GAMMa}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+a)}{(b + \xi_{ir})^{r+a}}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+a-1)}{(b + \xi_{ir})^{r+a-1}}}
	$
	
	\textbf{Entropy Loss Function (ELF)}
	
	The Bayes estimator under uniform prior is:
	$
	\hat{V}_{\text{ENT\_UNI}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+1)}{\xi_{ir}^{r+1}}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r)}{\xi_{ir}^{r}}}
	$
	
	The Bayes estimator under Jeffreys prior is:
	$
	\hat{V}_{\text{ENT\_Jeffreys}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r)}{\xi_{ir}^{r+1}}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r-1)}{\xi_{ir}^{r-1}}}
	$
	
	The Bayes estimator under exponential prior is:
	$
	\hat{V}_{\text{ENT\_EXP}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+1)}{(w + \xi_{ir})^{r+1}}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r)}{(w + \xi_{ir})^{r}}}
	$
	
	The Bayes estimator under gamma prior is:
	$
	\hat{V}_{\text{ENT\_Gamma}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+a)}{(b + \xi_{ir})^{r+a}}}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+a-1)}{(b + \xi_{ir})^{r+a-1}}}
	$
	
	Where $\psi(x) = \frac{\Gamma'(x)}{\Gamma(x)}$ is the digamma function.
	
	\textbf{Squared-Log Error Loss Function (SLELF)}
	
	The Bayes estimator under uniform prior is:
	$
	\hat{V}_{\text{SLELF\_UNI}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+1)}{\xi_{ir}^{r+1}} \left(\frac{\exp(\psi(r+1))}{\xi_{ir}}\right)}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+1)}{\xi_{ir}^{r+1}}}
	$
	
	The Bayes estimator under Jeffreys prior is:
	$
	\hat{V}_{\text{SLELF\_Jeffreys}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r)}{\xi_{ir}^{r}} \left(\frac{\exp(\psi(r))}{\xi_{ir}}\right)}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r)}{\xi_{ir}^{r}}}
	$
	
	The Bayes estimator under exponential prior is:
	$
	\hat{V}_{\text{SLELF\_EXP}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+1)}{(w + \xi_{ir})^{r+1}} \left(\frac{\exp(\psi(r+1))}{w + \xi_{ir}}\right)}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+1)}{(w + \xi_{ir})^{r+1}}}
	$
	
	The Bayes estimator under gamma prior is:
	$
	\hat{V}_{\text{SLELF\_gamma}} = \frac{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+a)}{(b + \xi_{ir})^{r+a}} \left(\frac{\exp(\psi(r+a))}{b + \xi_{ir}}\right)}{\sum_{k=0}^{n-r} (-1)^k \binom{n-r}{k} \frac{\Gamma(r+a)}{(b + \xi_{ir})^{r+a}}}
	$

## How to Use

The functions are implemented in R. To use them, you can simply copy the code into your R environment or script. Below are the functions implemented:

1. `PLF_loss_uniform`
2. `PLF_loss_Jeffreys`
3. `PLF_loss_exponential`
4. `PLF_loss_gamma`
5. `WLF_uniform`
6. `WLF_Jeffreys`
7. `WLF_exponential`
8. `WLF_gamma`
9. `ELF_uniform`
10. `ELF_Jeffreys`
11. `ELF_exponential`
12. `ELF_gamma`
13. `SLELF_uniform`
14. `SLELF_Jeffreys`
15. `SLELF_exponential`
16. `SLELF_gamma`

Each function takes parameters `x` (observations), `r`, `k`, and additional parameters specific to certain prior distributions.

## Example Usage

```R
# Load the functions
source("Bayesian_Loss_Functions.R")

# Example parameters
x <- c(1, 2, 3, 4, 5)
r <- 3
k <- 2
w <- 0.5
a <- 1
b <- 0.5

# Compute PLF estimators under different priors
V_PLF_U <- PLF_loss_uniform(x, r, k)
V_PLF_Jeffreys <- PLF_loss_Jeffreys(x, r, k)
V_PLF_EXP <- PLF_loss_exponential(x, r, k, w)
V_PLF_gamma <- PLF_loss_gamma(x, r, k, a, b)

# Print results
print(V_PLF_U)
print(V_PLF_Jeffreys)
print(V_PLF_EXP)
print(V_PLF_gamma)

# Similarly, compute estimators for other loss functions and priors
```

## Reference


These functions are based on Bayesian estimation principles and loss functions commonly used in statistical modeling. If you use these functions in your work, consider citing the original source or this repository.

**Reference**:
Sindhu, T. N., Saleem, M., & Aslam, M. (2013). Bayesian estimation for Topp Leone distribution under trimmed samples. Journal of Basic and Applied Scientific Research, 3(10), 347-360.

This publication provides insights into Bayesian estimation methods for the Topp Leone distribution, which serves as a basis for the implementation of the loss functions presented here.
