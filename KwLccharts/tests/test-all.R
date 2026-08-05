library(KwLccharts)

cat("Test 1: Lambert W inversion ... ")
stopifnot(abs(lambertW0(1) - 0.5671433) < 1e-6)
w <- lambertWm1(-0.2)
stopifnot(abs(w * exp(w) - (-0.2)) < 1e-8)
cat("OK\n")

cat("Test 2: KwL density integrates to 1 ... ")
val <- stats::integrate(function(x) dKwL(x, theta = 1.2, a = 1.5, b = 0.8), 0, Inf)$value
stopifnot(abs(val - 1) < 1e-4)
cat("OK\n")

cat("Test 3: KwL quantile inverts the cdf ... ")
us <- c(0.01, 0.1, 0.5, 0.9, 0.99)
xs <- qKwL(us, theta = 1.2, a = 1.5, b = 0.8)
back <- pKwL(xs, theta = 1.2, a = 1.5, b = 0.8)
stopifnot(max(abs(back - us)) < 1e-6)
cat("OK\n")

cat("Test 4: rKwL generates finite, non-negative variates ... ")
set.seed(1)
y <- rKwL(500, theta = 1, a = 1.5, b = 1.5)
stopifnot(all(is.finite(y)), all(y >= 0))
cat("OK\n")

cat("Test 5: fit_KwL recovers parameters with large samples ... ")
set.seed(42)
y <- rKwL(20000, theta = 1, a = 1.5, b = 2)
fit <- fit_KwL(y)
stopifnot(abs(fit$theta - 1)   < 0.15)
stopifnot(abs(fit$a     - 1.5) < 0.15)
stopifnot(abs(fit$b     - 2)   < 0.3)
cat("OK\n")

cat("\nAll KwLccharts tests passed.\n")
