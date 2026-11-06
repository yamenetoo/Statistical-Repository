\name{bkwg_mixture}
\alias{bkwg_mixture}
\alias{dbkwg_series}
\title{Infinite-Mixture Representation of the BKw-G Density}
\description{
The BKw-G density is a linear combination of Kumaraswamy-G densities:
\deqn{f(x) = \sum_{j=0}^\infty w_j \, f_{KwG}(x; a, b(\beta + j)),}{f(x) = sum_j w_j f_KwG(x; a, b (beta + j)),}
with
\deqn{w_j = \frac{(-1)^j \binom{\alpha - 1}{j}}{B(\alpha, \beta)(\beta + j)}.}{w_j = (-1)^j choose(alpha - 1, j) / (B(alpha, beta) (beta + j)).}
\code{bkwg_mixture} returns the weights and \code{dbkwg_series} evaluates the
truncated series.
}
\usage{
bkwg_mixture(alpha, beta, J = 100)
dbkwg_series(x, a = 1, b = 1, alpha = 1, beta = 1,
             baseline = "weibull", par = NULL, J = 100)
}
\arguments{
  \item{alpha, beta}{positive beta-generator shape parameters.}
  \item{J}{truncation point of the series (the series is finite when
    \code{alpha} is a positive integer).}
  \item{x}{vector of quantiles.}
  \item{a, b}{positive Kumaraswamy shape parameters.}
  \item{baseline, par}{baseline distribution and its parameters, as in
    \code{\link{dbkwg}}.}
}
\details{
The representation follows from the generalised binomial expansion of
\eqn{K(x)^{\alpha - 1}}{K(x)^(alpha - 1)}. It allows moments and other
properties of the BKw-G family to be obtained from those of the
Kumaraswamy-G family. The weights sum to one.
}
\value{
\code{bkwg_mixture} returns a data frame with columns \code{j},
\code{weight} and \code{b_multiplier} (\eqn{\beta + j}{beta + j}; the
\eqn{j}-th component is Kw-G with parameters \eqn{a} and
\eqn{b(\beta + j)}{b (beta + j)}). \code{dbkwg_series} returns the
approximated density.
}
\seealso{\code{\link{dbkwg}}}
\examples{
w <- bkwg_mixture(alpha = 2.5, beta = 1.5, J = 200)
sum(w$weight)
x <- c(0.5, 1, 2)
dbkwg_series(x, a = 2, b = 1, alpha = 2.5, beta = 1.5, par = c(2, 1))
dbkwg(x, a = 2, b = 1, alpha = 2.5, beta = 1.5, par = c(2, 1))
}
\keyword{distribution}
