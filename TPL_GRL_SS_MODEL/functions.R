#tpl cdf
F_x <- function(x, v) {  
  # Nadarajah & S. Kotz
  if (all(x > 0 & x <= 1) && all(v > 0 & v < 1)) {
    ptopple(x, v, lower.tail = TRUE, log.p = FALSE)
  } else {
    stop("Values of x and v must be between 0 and 1")
  }
}
#GRL pdf
g_Y <- function(x, alpha, beta) {
  if (all(alpha > 0) && all(beta > 0) && all(x > 0 )) {
    dgenray(x, scale = 1/alpha, shape = beta, log = FALSE)
  } else {
    stop("Values of alpha and beta must be greater than 0, and values of x must be between 0 and 1")
  }
}


integrand <- function(x, alpha, beta, v,c_value=1) {
  (1 - F_x(x, v)) * g_Y(c_value*x, alpha, beta)
}
# R 
R <- function(theta1, theta2,c_value=1) {
  v <- theta1[1]
  alpha <- theta2[1]
  beta <- theta2[2]
  # Perform numerical integration
  result <- integral(integrand, 1E-5, 1, alpha = alpha, beta = beta, v = v,,c_value)
  return(result)
}

############################################################
######Reliability in Multi-Stress Environments##############
############################################################
Multi_Stress_Environments <- function(k, Theta1, Theta2) {
  results <- 1
  for (i in 1:k) {
    results <-results* R(Theta1, as.numeric(Theta2[i,]))
  }
  return(results)
}
#Example
#Theta2=matrix(c(0.1,0.1,0.2,0.2),2,2)
#Multi_Stress_Environments(2,0.1,Theta2)

############################################################
##########################R3################################
############################################################

R3 <- function(theta1, theta2, theta3, c_value = 1){
  R_theta1_theta3 <- R(theta1, theta3, c_value)
  R_theta2_theta3 <- R(theta2, theta3, c_value)
  result <- (1 - R_theta2_theta3) * R_theta1_theta3
  return(result)
}
#X~theta1 
#Z~theta2 
#Y~theta3 
#R3(0.7,0.3,c(0.1,0.1))


###########################
########## casecade #######
###########################
Ri_casecade <- function(i, Theta1, Theta2, c_value) {
  if (i == 1) {
    theta1_1 <- as.numeric(Theta1[1])
    theta2_1 <- as.numeric(Theta2[1, ])
    return(R(theta1_1, theta2_1))
  }
  
  # Define Theta2_star values
  Theta2_star_alpha <- Theta2[i, "alpha"]/ c_value^(i-1)
  Theta2_star_beta <- Theta2[i, "beta"] 
  
  R_values <- 1
  for (j in 1:(i-1))
  {
    theta1_j <- as.numeric(Theta1[j])
    theta2_j_alpha <- Theta2[j, "alpha"]/ c_value^(j-1)
    theta2_j_beta <- Theta2[j, "beta"] 
    
    R_values <- R_values * (1 - R(theta1_j, c(theta2_j_alpha, theta2_j_beta),c_value^(j-1)))
  }
  theta1_i <- as.numeric(Theta1[i])
  R_values <- R_values * R(theta1_i, c(Theta2_star_alpha, Theta2_star_beta),c_value^(i-1))
  
  return(R_values)
}

casecade=function(n, Theta1, Theta2,c_value){
  Ri_values <- c()
  for (i in 1:n) {
    Ri_values[i] <- Ri_casecade(i, Theta1, Theta2,c_value)
  }
  return(sum(Ri_values))
}


