

MC=1000
N <- c(10, 25, 50)
Theta1 <- c(0.3, 0.5, 0.7, 0.9)
Theta2 <- matrix(c(
  0.5, 0.5,
  1, 1,
  0.5, 1,
  2, 2,
  3, 3.5,
  4.5, 1
), ncol = 2, byrow = TRUE)

TPL_mle <-  numeric(length(N) * length(Theta1) * nrow(Theta2)*MC)
TPL_JacK <- numeric(length(N) * length(Theta1) * nrow(Theta2)*MC)
TPL_McMc<-  numeric(length(N) * length(Theta1) * nrow(Theta2)*MC)
GRL_Jack1 <-numeric(length(N) * length(Theta1) * nrow(Theta2)*MC)
GRL_Jack2 <-numeric(length(N) * length(Theta1) * nrow(Theta2)*MC)
GRL_MLE1 <- numeric(length(N) * length(Theta1) * nrow(Theta2)*MC)
GRL_MLE2 <- numeric(length(N) * length(Theta1) * nrow(Theta2)*MC)
GRL_MCMC1 <-numeric(length(N) * length(Theta1) * nrow(Theta2)*MC)
GRL_MCMC2 <-numeric(length(N) * length(Theta1) * nrow(Theta2)*MC)

par <- character(length(N) * length(Theta1) * nrow(Theta2)*MC)
i <- 0
for (n in N) {
  for (theta1 in Theta1) {
    for (j in 1:nrow(Theta2)) {
        for (rep in 1:MC){
        i <- i + 1
        theta2 <- Theta2[j,]
        theta21 <- Theta2[j,1]
        theta22 <- Theta2[j,2]
        par[i] <- paste(n, theta1, theta21, theta22)
        tryCatch({
          Y <- rgenray(n = n, scale = theta2[1], shape = theta2[2])
          X <- rtopple(n = n, shape = theta1)
          
          TPL_mle [i] <- TplMle(X)
          TPL_JacK[i] <- TplJack(X)
          TPL_McMc[i] <- TPL_Baysian(X)
          
          res          <- GrlJack(Y)
          GRL_Jack1[i] <- res[1]
          GRL_Jack2[i] <- res[2]
          
          res         <- GrlMle(Y)
          GRL_MLE1[i] <- res[1]
          GRL_MLE2[i] <- res[2]
          res         <- GRL_Baysian(Y)
          GRL_MCMC1[i]<- res[1]
          GRL_MCMC2[i]<- res[2]
    
          print(paste(par[i],TPL_McMc[i]))
        },error = function(e){
          cat("An error occurred:", conditionMessage(e), "\n")
        })
      }
    }
  }
}
################################################################################
################################################################################
file_path <- "combined_matrix.csv"
par_split <- strsplit   (par,"  ")
par_matrix <- do.call   (rbind, par_split)
par_matrix <- apply     (par_matrix, 2, as.numeric)
colnames(par_matrix)=c  ("n","v","alpha","beta")
combined_matrix <- cbind(par      ,par_matrix, TPL_mle   ,
                         TPL_JacK ,TPL_McMc  , GRL_Jack1 ,
                         GRL_Jack2,GRL_MLE1  , GRL_MLE2  ,
                         GRL_MCMC1,GRL_MCMC2)
write.csv(combined_matrix, file = file_path, row.names = FALSE)
cat("CSV file saved:", file_path, "\n")
