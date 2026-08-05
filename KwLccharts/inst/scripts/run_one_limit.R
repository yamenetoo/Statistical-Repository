args <- commandArgs(trailingOnly = TRUE)
# args: n u theta a b B MC outfile
n     <- as.integer(args[1])
u     <- as.numeric(args[2])
theta <- as.numeric(args[3])
a     <- as.numeric(args[4])
b     <- as.numeric(args[5])
B     <- as.integer(args[6])
MC    <- as.integer(args[7])
outfile <- args[8]

library(KwLccharts)
set.seed(1000 + n * 100 + round(u * 1000) + round(theta*10) + round(a*10) + round(b*10))

alpha <- 0.0027
t0 <- Sys.time()
res <- KwLBootCC(u = u, alpha = alpha, theta = theta, a = a, b = b, n = n, B = B, MC = MC)
t1 <- Sys.time()

row <- data.frame(n = n, theta = theta, a = a, b = b, u = u,
                   MLCL = res$MLCL, SDLCL = res$SDLCL, MUCL = res$MUCL, SDUCL = res$SDUCL,
                   B = B, MC = MC)

if (!file.exists(outfile)) {
  write.csv(row, outfile, row.names = FALSE)
} else {
  write.table(row, outfile, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
}

cat(sprintf("n=%d u=%.2f theta=%.2f a=%.2f b=%.2f -> MLCL=%.5f SDLCL=%.5f MUCL=%.5f SDUCL=%.5f  [%.1fs]\n",
            n, u, theta, a, b, res$MLCL, res$SDLCL, res$MUCL, res$SDUCL, as.numeric(t1 - t0)))
