results_df <- data.frame()
n_values <- c(4,5,6)
u_values <- c(0.01, 0.1, 0.9, 0.99)
p_values <- c(0.1, 0.5, 0.9)
theta_values <- c(1, 2)
for (n in n_values)
{
  for (u in u_values)
  {
    for (p in p_values)
    {
      for (theta in theta_values)
      {
        result <- LGBootStrapControlChart(u = u, alpha = 0.07, theta = theta, p = p, n = n,
                                          BootStrap_iterations = 10000, MontiCarlo_iter = 10,
                                          ARL_rep = 1000, verbos = TRUE)
        cat("n=",n,"p",p,"th=",theta,"\n")
        result_df <- data.frame(u = u, p = p, theta = theta,
                                MLCL = result$MLCL, MUCL = result$MUCL,
                                SDLCL = result$SDLCL, SDUCL = result$SDUCL,
                                MARL = result$MARL, SDARL = result$SDARL)
        
        results_df <- rbind(results_df, result_df)
      }
    }
  }
}


write.csv(results_df, "res.csv", row.names = FALSE)
