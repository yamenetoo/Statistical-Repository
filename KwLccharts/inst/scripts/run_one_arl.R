args <- commandArgs(trailingOnly = TRUE)
# args: u LCL UCL theta1 a1 b1 label n_runs max_run outfile
u       <- as.numeric(args[1])
LCL     <- as.numeric(args[2])
UCL     <- as.numeric(args[3])
theta1  <- as.numeric(args[4])
a1      <- as.numeric(args[5])
b1      <- as.numeric(args[6])
label   <- args[7]
n_runs  <- as.integer(args[8])
max_run <- as.integer(args[9])
outfile <- args[10]

library(KwLccharts)
set.seed(2000 + round(u*1000) + round(theta1*10) + round(a1*10) + round(b1*10))

t0 <- Sys.time()
res <- KwL_ARL(u, LCL, UCL, theta1, a1, b1, n = 5, n_runs = n_runs, max_run = max_run)
t1 <- Sys.time()

row <- data.frame(u = u, shift = label, theta1 = theta1, a1 = a1, b1 = b1,
                   ARL = res$ARL, SDRL = res$SDRL, n_runs = n_runs, max_run = max_run)

if (!file.exists(outfile)) {
  write.csv(row, outfile, row.names = FALSE)
} else {
  write.table(row, outfile, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
}

cat(sprintf("u=%.2f [%s] theta1=%.2f a1=%.2f b1=%.2f -> ARL=%.3f SDRL=%.3f  [%.1fs]\n",
            u, label, theta1, a1, b1, res$ARL, res$SDRL, as.numeric(difftime(t1, t0, units = "secs"))))
