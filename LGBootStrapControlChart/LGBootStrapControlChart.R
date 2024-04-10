if (!requireNamespace("LindleyPowerSeries", quietly = TRUE))
  install.packages("LindleyPowerSeries")

if (!requireNamespace("fitdistrplus", quietly = TRUE))
  install.packages("fitdistrplus")

library(LindleyPowerSeries)
library(fitdistrplus)
options(show.error.messages = FALSE)
rm(list = ls())


dlindleygeometric <- function(x, theta, p)
{
  if (any(is.na(x) | x <= 0 | theta <= 0 | p < 0 | p > 1))
  {
    return(rep(NaN, length(x)))
  }
  return(LindleyPowerSeries::dlindleygeometric(x, theta, p))
}

plindleygeometric <- function(q, theta, p)
{
  if (any(is.na(q) | q <= 0 | theta <= 0 | p < 0 | p > 1)) {
    return(rep(NaN, length(q)))
  }
  return(LindleyPowerSeries::plindleygeometric(q, theta, p))
}


fit_distribution <- function(sample)
{
  result <- tryCatch(
    {
      fitlingeo <- fitdist(data = sample,
                           distr = "lindleygeometric",
                           start = list(theta = 0.5, p = 0.01),
                           method = "mle",
                           control = list(trace = FALSE),
                           optim.method = "Nelder-Mead")
      para <- as.vector(summary(fitlingeo)$estimate)
      phat <- para[1]
      thetahat <- para[2]
      #cat("Estimated parameters:", paste(phat, thetahat), "\n")
      return(c(phat, thetahat))
    },
    error = function(e)
    {
      return(NULL)
    })
  return(result)
}



LGBootStrapControlChart <- function(u, alpha, theta, p, n=5,
                                    BootStrap_iterations=1E+2,
                                    MontiCarlo_iter=30,
                                    ARL_rep = 100,verbos=FALSE)
{
  QS<-LCL <- UCL <- ARL_ <- NULL
  
  for(itr in 1:MontiCarlo_iter)
  {
    P<- Theta <-NULL
    iteration<-1
    while (iteration <= BootStrap_iterations)
    {
      sample <- rlindleygeometric(n, theta, p)
      result <- fit_distribution(sample)
      if (!is.null(result))
      {
        Theta[iteration]=result[1]
        P[iteration]=result[2]
        iteration <- iteration + 1
      }
    }
    for (i in 1:BootStrap_iterations)   QS[i] <- qlindleygeometric(u, Theta[i], P[i])
    QS <- QS[order(QS)]
    LCL[itr] <- QS[as.integer(BootStrap_iterations * alpha / 2)]
    UCL[itr] <- QS[as.integer(BootStrap_iterations * (1 - alpha / 2))]
    if(verbos)        cat("LCL[",itr,"]=", LCL[itr], "\tUCL[",itr,"]=", UCL[itr], "\n")
  }
  
  MLCL <- mean(LCL)
  MUCL <- mean(UCL)
  SDLCL <- var(LCL)**0.5
  SDUCL <- var(UCL)**0.5
  if(verbos)
  {
    cat("###############################################################################################\n")
    cat("Mean of LCL (MLCL): ", MLCL,"\t  Mean of UCL (MUCL): ", MUCL, "\n")
    cat("Standard deviation of LCL (SDLCL): ", SDLCL,"\t  Standard deviation of UCL (SDUCL): ", SDUCL, "\n")
    cat("###############################################################################################\n")
  }
  
  for(itr in 1:ARL_rep)
  {
    ARL <- 1
    continue <- TRUE
    while(continue)
    {
      sample <- rlindleygeometric(n, theta, p)
      result <- fit_distribution(sample)
      if (!is.null(result))
      {
        theta_hat <- result[1]
        p_hat <- result[2]
        Q_hat <- qlindleygeometric(u, theta_hat, p_hat)
        if (Q_hat < MLCL || Q_hat > MUCL)
        {
          ARL_[itr] <- ARL
          continue <- FALSE
          #if(verbos)  cat("ARL[", itr, "]=", ARL, "\n")
        }else
        {
          ARL <- ARL + 1
        }
      }
    }
  }
  if(verbos)
  {
    cat("#############################################################################################\n")
    MARL=mean(ARL_)       ;cat("Mean of ARL (MARL): ", MARL , "\n")
    SDARL=var(ARL_)**0.5  ;cat("SD of ARL (SDARL): " , SDARL, "\n")
    cat("###############################################################################################\n")
  }
  return(list(MLCL=MLCL,MUCL=MUCL,SDLCL=SDLCL,SDUCL=SDUCL,MARL=MARL,SDARL=SDARL))
}
