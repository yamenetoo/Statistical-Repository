\name{bkwg_order}
\alias{bkwg_order}
\alias{dos_bkwg}
\alias{pos_bkwg}
\title{Order Statistics of the BKw-G Distribution}
\description{
Density and distribution function of the \eqn{i}-th order statistic in a
random sample of size \eqn{N} from a BKw-G distribution.
}
\usage{
dos_bkwg(x, i, N, a = 1, b = 1, alpha = 1, beta = 1,
         baseline = "weibull", par = NULL, log = FALSE)
pos_bkwg(q, i, N, a = 1, b = 1, alpha = 1, beta = 1,
         baseline = "weibull", par = NULL)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{i}{order of the statistic, \eqn{1 \le i \le N}{1 <= i <= N}.}
  \item{N}{sample size.}
  \item{a, b, alpha, beta}{positive shape parameters.}
  \item{baseline, par}{baseline distribution and its parameters, as in
    \code{\link{dbkwg}}.}
  \item{log}{logical; return the log-density.}
}
\details{
\deqn{f_{i:N}(x) = \frac{N!}{(i-1)!(N-i)!} f(x) F(x)^{i-1} [1 - F(x)]^{N-i},}{f_i:N(x) = N! / ((i-1)! (N-i)!) f(x) F(x)^(i-1) [1 - F(x)]^(N-i),}
and \eqn{F_{i:N}(x) = I_{F(x)}(i, N - i + 1)}{F_i:N(x) = I_F(x)(i, N - i + 1)}.
}
\value{Numeric vector.}
\seealso{\code{\link{dbkwg}}}
\examples{
## density of the sample minimum and maximum for N = 10
x <- seq(0.05, 3, by = 0.05)
fmin <- dos_bkwg(x, 1, 10, a = 2, b = 1, alpha = 1.5, beta = 2, par = c(2, 1))
fmax <- dos_bkwg(x, 10, 10, a = 2, b = 1, alpha = 1.5, beta = 2, par = c(2, 1))
pos_bkwg(1, 5, 10, a = 2, b = 1, alpha = 1.5, beta = 2, par = c(2, 1))
}
\keyword{distribution}
