# Maximum Likelihood Estimation (MLE) Functions
TplMle <- function(data_) {
  tdata <- data.frame(y = data_)
  tfit <- vglm(y ~ 1, topple, tdata, trace = FALSE, crit = "coef")
  return(as.numeric(Coef(tfit)))
}

GrlMle <- function(data_) {
  data_ <- data.frame(y = data_)
  fit <- vglm(y ~ 1, genrayleigh, data = data_, trace = FALSE)
  return(as.numeric(Coef(fit)))
}

# Jackknife Estimation Functions
TplJack <- function(data) {
  n <- length(data)
  statistics <- numeric(n)
  for (i in 1:n) {
    statistics[i] <- TplMle(data[-i])  # Leave out observation i
  }
  res <- mean(statistics)
  return(res)
}

GrlJack <- function(y) {
  n <- length(y)  
  statistics <- matrix(0, as.numeric(n), 2)
  for (i in 1:n) {
    statistics[i,] <- GrlMle(y[-i])  # Leave out column i
  }
  
  res <- colMeans(statistics)
  return(res)
}

# MCMC Functions
TPL_Baysian <- function(x, initial_v = 0.5, shape = 1, rate = 1, n_iterations = 10000, n_burnin = 1000) {
  likelihood_tpl <- function(v, x) {
    if (v <= 0) return(rep(0, length(x)))
    invalid_indices <- x <= 0 | x >= 2
    if (any(invalid_indices)) return(rep(0, length(x)))
    likelihood_values <- 2 * v * x^(v-1) * (1 - x) * (2 - x)^(v-1)
    likelihood_values[invalid_indices] <- 0
    return(likelihood_values)
  }
  prior_v <- function(v, shape, rate) {
    return(dgamma(v, shape = shape, rate = rate))
  }
  
  log_posterior_tpl <- function(params, x, shape, rate) {
    v <- params
    if (v <= 0) return(-Inf)
    
    prior_prob_v <- log(prior_v(v, shape, rate))
    if (!is.finite(prior_prob_v)) return(-Inf)
    
    likelihood_values <- likelihood_tpl(v, x)
    if (any(!is.finite(likelihood_values))) return(-Inf)
    
    invalid_indices <- which(likelihood_values == 0)
    if (length(invalid_indices) > 0) return(-Inf)
    
    log_posterior <- sum(prior_prob_v) + sum(log(likelihood_values))
    if (!is.finite(log_posterior)) return(-Inf)
    
    return(log_posterior)
  }
  
  mcmc_result_tpl <- MCMCmetrop1R(log_posterior_tpl, 
                                  theta.init = initial_v, 
                                  x = x, 
                                  shape = shape, 
                                  rate = rate, 
                                  burnin = n_burnin, 
                                  thin = 1,
                                  verbose = FALSE
  )
  res <- summary(mcmc_result_tpl)
  res <- res$statistics
  v_mcmc <- res[1]
  if (v_mcmc > 1) v_mcmc <- 0.999
  return(v_mcmc)
}

GRL_Baysian <- function(y, shape_alpha = 1, rate_alpha = 1, shape_lambda = 1, rate_lambda = 1, n_iterations = 10000, n_burnin = 1000) {
  likelihood <- function(alpha, lambda, y) {
    if (alpha <= 0 | lambda <= 0) return(rep(0, length(y)))
    invalid_indices <- y <= 0
    if (any(invalid_indices)) return(rep(0, length(y)))
    likelihood_values <- 2 * alpha * lambda^2 * y * exp(-(lambda * y)^2) * (1 - exp(-(lambda * y)^2))^(alpha - 1)
    likelihood_values[invalid_indices] <- 0
    return(likelihood_values)
  }
  
  prior_alpha <- function(alpha, shape_alpha, rate_alpha) {
    return(dgamma(alpha, shape = shape_alpha, rate = rate_alpha))
  }
  
  prior_lambda <- function(lambda, shape_lambda, rate_lambda) {
    return(dgamma(lambda, shape = shape_lambda, rate = rate_lambda))
  }
  
  log_posterior <- function(params, y, shape_alpha, rate_alpha, shape_lambda, rate_lambda) {
    alpha <- params[1]
    lambda <- params[2]
    
    if (alpha <= 0 | lambda <= 0) return(-Inf)
    
    prior_prob_alpha <- log(prior_alpha(alpha, shape_alpha, rate_alpha))
    prior_prob_lambda <- log(prior_lambda(lambda, shape_lambda, rate_lambda))
    
    if (!is.finite(prior_prob_alpha) | !is.finite(prior_prob_lambda)) return(-Inf)
    
    likelihood_values <- likelihood(alpha, lambda, y)
    if (any(!is.finite(likelihood_values))) return(-Inf)
    
    invalid_indices <- which(likelihood_values == 0)
    if (length(invalid_indices) > 0) return(-Inf)
    
    log_posterior <- sum(prior_prob_alpha) + sum(prior_prob_lambda) + sum(log(likelihood_values))
    if (!is.finite(log_posterior)) return(-Inf)
    
    return(log_posterior)
  }
  
  initial_values <- c(alpha = 1, lambda = 1)
  
  mcmc_result <- MCMCmetrop1R(log_posterior,
                              theta.init = initial_values, 
                              y = y,
                              shape_alpha = shape_alpha, 
                              rate_alpha = rate_alpha, 
                              shape_lambda = shape_lambda, 
                              rate_lambda = rate_lambda, 
                              burnin = n_burnin, 
                              thin = 1,
                              verbose = FALSE
  )
  res <- summary(mcmc_result)
  res <- res$statistics[,1]
  res <- c(1/res[2], res[1])
  return(res)
}
