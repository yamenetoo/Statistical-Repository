
if (!requireNamespace("gtools", quietly = TRUE)) {
  install.packages("gtools")
}

if (!requireNamespace("gsl", quietly = TRUE)) {
  install.packages("gsl")
}


if (!requireNamespace("VGAM", quietly = TRUE)) {
  install.packages("VGAM")
}




library(gtools)
library(gsl) 
library(VGAM)



xi <- function(x, r, k)
{

  # Sort the input vector x in ascending order
  x_sorted <- sort(x)
  
  # Compute the sum of logs for the first r values and the r-th value
  sum_logs <- sum(log(2 * x_sorted[1:r] - x_sorted[1:r]^2))
  rth_log <- k * log(2 * x_sorted[r] - x_sorted[r]^2)
  
  # Compute xi_ir
  xi_ir <- -(sum_logs + rth_log)
  
  return(xi_ir)
}


###########################
# PLF loss function under #
###########################

PLF_loss_uniform <- function(x, r, k) {
  xi_ir <- xi(x, r, k)
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 3) / (xi_ir^(r + 3))
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 1) / (xi_ir^(r + 1))
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_PLF_U_hat <- sqrt(numerator / denominator)
  
  return(V_PLF_U_hat)
}

# PLF loss function under Jeffreys prior
PLF_loss_Jeffreys <- function(x, r, k) {
  xi_ir <- xi(x, r, k)
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 2) / (xi_ir^(r + 2))
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r) / (xi_ir^r)
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_PLF_Jeffreys_hat <- sqrt(numerator / denominator)
  
  return(V_PLF_Jeffreys_hat)
}

# PLF loss function under exponential prior
PLF_loss_exponential <- function(x, r, k, w) {
  xi_ir <- xi(x, r, k)
  xi_ir_w <- xi_ir + w
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 3) / (xi_ir_w^(r + 3))
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 1) / (xi_ir_w^(r + 1))
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_PLF_EXP_hat <- sqrt(numerator / denominator)
  
  return(V_PLF_EXP_hat)
}

# PLF loss function under gamma prior
PLF_loss_gamma <- function(x, r, k, a, b) {
  xi_ir <- xi(x, r, k)
  xi_ir_b <- xi_ir + b
  a_r <- a + r
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(a_r + 2) / (xi_ir_b^(a_r + 2))
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(a_r) / (xi_ir_b^a_r)
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_PLF_gamma_hat <- sqrt(numerator / denominator)
  
  return(V_PLF_gamma_hat)
}


################################
# Weighted Loss Function (WLF) #
################################

WLF_uniform <- function(x, r, k) {
  xi_ir <- xi(x, r, k)
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 1) / (xi_ir^(r + 1))
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r) / (xi_ir^r)
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_WLF_UNI_hat <- numerator / denominator
  
  return(V_WLF_UNI_hat)
}

# Weighted Loss Function (WLF) under Jeffreys prior
WLF_Jeffreys <- function(x, r, k) {
  xi_ir <- xi(x, r, k)
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r) / (xi_ir^r)
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r - 1) / (xi_ir^(r - 1))
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_WLF_Jeffreys_hat <- numerator / denominator
  
  return(V_WLF_Jeffreys_hat)
}

# Weighted Loss Function (WLF) under exponential prior
WLF_exponential <- function(x, r, k, w) {
  xi_ir <- xi(x, r, k)
  xi_ir_w <- xi_ir + w
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 1) / (xi_ir_w^(r + 1))
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r) / (xi_ir_w^r)
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_WLF_EXP_hat <- numerator / denominator
  
  return(V_WLF_EXP_hat)
}

# Weighted Loss Function (WLF) under gamma prior
WLF_gamma <- function(x, r, k, a, b) {
  xi_ir <- xi(x, r, k)
  xi_ir_b <- xi_ir + b
  a_r <- a + r
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(a_r) / (xi_ir_b^(a_r))
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(a_r - 1) / (xi_ir_b^(a_r - 1))
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_WLF_Gamma_hat <- numerator / denominator
  
  return(V_WLF_Gamma_hat)
}





###############################
# Entropy Loss Function (ELF) #
###############################
#Entropy Loss Function (ELF) under uniform prior
ELF_uniform <- function(x, r, k) {
  xi_ir <- xi(x, r, k)
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 1) / (xi_ir^(r + 1))
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r) / (xi_ir^r)
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_ENT_UNI_hat <- numerator / denominator
  
  return(V_ENT_UNI_hat)
}

# Entropy Loss Function (ELF) under Jeffreys prior
ELF_Jeffreys <- function(x, r, k) {
  xi_ir <- xi(x, r, k)
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r) / (xi_ir^(r + 1))
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r - 1) / (xi_ir^(r - 1))
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_ENT_Jeffreys_hat <- numerator / denominator
  
  return(V_ENT_Jeffreys_hat)
}

# Entropy Loss Function (ELF) under exponential prior
ELF_exponential <- function(x, r, k, w) {
  xi_ir <- xi(x, r, k)
  xi_ir_w <- xi_ir + w
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 1) / (xi_ir_w^(r + 1))
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r) / (xi_ir_w^r)
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_ENT_EXP_hat <- numerator / denominator
  
  return(V_ENT_EXP_hat)
}

# Entropy Loss Function (ELF) under gamma prior
ELF_gamma <- function(x, r, k, a, b) {
  xi_ir <- xi(x, r, k)
  xi_ir_b <- xi_ir + b
  a_r <- a + r
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(a_r) / (xi_ir_b^(a_r))
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(a_r - 1) / (xi_ir_b^(a_r - 1))
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_ENT_Gamma_hat <- numerator / denominator
  
  return(V_ENT_Gamma_hat)
}



############################################
## Squared-Log Error Loss Function (SLELF)##
############################################

# Squared-Log Error Loss Function (SLELF) under uniform prior
SLELF_uniform <- function(x, r, k) {
  xi_ir <- xi(x, r, k)
  psi_r_plus_1 <- digamma(r + 1)
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 1) / (xi_ir^(r + 1)) * exp(psi_r_plus_1) / xi_ir
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 1) / (xi_ir^(r + 1))
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_SLELF_UNI_hat <- numerator / denominator
  
  return(V_SLELF_UNI_hat)
}

# Squared-Log Error Loss Function (SLELF) under Jeffreys prior
SLELF_Jeffreys <- function(x, r, k) {
  xi_ir <- xi(x, r, k)
  psi_r <- digamma(r)
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r) / (xi_ir^r) * exp(psi_r) / xi_ir
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r) / (xi_ir^r)
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_SLELF_Jeffreys_hat <- numerator / denominator
  
  return(V_SLELF_Jeffreys_hat)
}

# Squared-Log Error Loss Function (SLELF) under exponential prior
SLELF_exponential <- function(x, r, k, w) {
  xi_ir <- xi(x, r, k)
  xi_ir_w <- xi_ir + w
  psi_r_plus_1 <- digamma(r + 1)
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 1) / (xi_ir_w^(r + 1)) * exp(psi_r_plus_1) / xi_ir_w
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(r + 1) / (xi_ir_w^(r + 1))
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_SLELF_EXP_hat <- numerator / denominator
  
  return(V_SLELF_EXP_hat)
}

# Squared-Log Error Loss Function (SLELF) under gamma prior
SLELF_gamma <- function(x, r, k, a, b) {
  xi_ir <- xi(x, r, k)
  xi_ir_b <- xi_ir + b
  a_r <- a + r
  psi_a_r <- digamma(a_r)
  
  numerator <- 0
  denominator <- 0
  
  for (k in 0:(length(x) - r)) {
    numerator_term <- (-1)^k * choose(length(x) - r, k) * gamma(a_r) / (xi_ir_b^(a_r)) * exp(psi_a_r) / xi_ir_b
    denominator_term <- (-1)^k * choose(length(x) - r, k) * gamma(a_r) / (xi_ir_b^(a_r))
    
    numerator <- numerator + numerator_term
    denominator <- denominator + denominator_term
  }
  
  V_SLELF_Gamma_hat <- numerator / denominator
  
  return(V_SLELF_Gamma_hat)
}


########################################################################################
########################################################################################
########################################################################################

# Examples:
#x <- 






r <- 3  # Example value of r
k <- 2  # Example value of k
w <- 0.5  # Example value of w
a <- 1  # Example value of a
b <- 0.5  # Example value of b

xi_value <- xi(x, r, k)
print(xi_value)


V_PLF_U <- PLF_loss_uniform(x, r, k)
V_PLF_Jeffreys <- PLF_loss_Jeffreys(x, r, k)
V_PLF_EXP <- PLF_loss_exponential(x, r, k, w)
V_PLF_gamma <- PLF_loss_gamma(x, r, k, a, b)

print(V_PLF_U)
print(V_PLF_Jeffreys)
print(V_PLF_EXP)
print(V_PLF_gamma)



V_WLF_UNI <- WLF_uniform(x, r, k)
V_WLF_Jeffreys <- WLF_Jeffreys(x, r, k)
V_WLF_EXP <- WLF_exponential(x, r, k, w)
V_WLF_Gamma <- WLF_gamma(x, r, k, a, b)

print(V_WLF_UNI)
print(V_WLF_Jeffreys)
print(V_WLF_EXP)
print(V_WLF_Gamma)


V_ENT_UNI <- ELF_uniform(x, r, k)
V_ENT_Jeffreys <- ELF_Jeffreys(x, r, k)
V_ENT_EXP <- ELF_exponential(x, r, k, w)
V_ENT_Gamma <- ELF_gamma(x, r, k, a, b)

print(V_ENT_UNI)
print(V_ENT_Jeffreys)
print(V_ENT_EXP)
print(V_ENT_Gamma)



V_SLELF_UNI <- SLELF_uniform(x, r, k)
V_SLELF_Jeffreys <- SLELF_Jeffreys(x, r, k)
V_SLELF_EXP <- SLELF_exponential(x, r, k, w)
V_SLELF_Gamma <- SLELF_gamma(x, r, k, a, b)

print(V_SLELF_UNI)
print(V_SLELF_Jeffreys)
print(V_SLELF_EXP)
print(V_SLELF_Gamma)
