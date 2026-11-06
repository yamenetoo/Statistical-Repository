\name{BKwG-package}
\alias{BKwG-package}
\alias{BKwG}
\docType{package}
\title{The Beta Kumaraswamy-G Family of Distributions}
\description{
Tools for the beta Kumaraswamy-G (BKw-G) family of distributions of
Handique, Chakraborty and Ali (2017). For a baseline distribution with cdf
\eqn{G(x)} and pdf \eqn{g(x)}, the BKw-G distribution has cdf
\deqn{F(x) = I_{K(x)}(\alpha, \beta), \qquad K(x) = 1 - [1 - G(x)^a]^b,}{F(x) = I_K(x)(alpha, beta),  K(x) = 1 - [1 - G(x)^a]^b,}
and pdf
\deqn{f(x) = \frac{a\, b\, g(x)\, G(x)^{a-1} [1 - G(x)^a]^{b\beta - 1} K(x)^{\alpha - 1}}{B(\alpha, \beta)},}{f(x) = a b g(x) G(x)^(a-1) [1 - G(x)^a]^(b beta - 1) K(x)^(alpha - 1) / B(alpha, beta),}
where \eqn{I_y(\alpha, \beta)}{I_y(alpha, beta)} is the regularised incomplete
beta function and \eqn{a, b, \alpha, \beta > 0}{a, b, alpha, beta > 0}.
The parameters \code{alpha} and \code{beta} correspond to \eqn{m} and
\eqn{n} in the notation of Handique et al. (2017).
}
\details{
Main functions:
\itemize{
  \item Distribution: \code{\link{dbkwg}}, \code{\link{pbkwg}},
    \code{\link{qbkwg}}, \code{\link{rbkwg}}, \code{\link{sbkwg}},
    \code{\link{hbkwg}}, \code{\link{chbkwg}}, \code{\link{score_bkwg}}.
  \item Baselines: \code{\link{bkwg_baselines}}, \code{\link{bkwg_baseline}}.
  \item Properties: \code{\link{mbkwg}}, \code{\link{mgf_bkwg}},
    \code{\link{bkwg_stats}}, \code{\link{entropy_bkwg}},
    \code{\link{mrl_bkwg}}, \code{\link{dos_bkwg}}, \code{\link{pos_bkwg}},
    \code{\link{bkwg_mixture}}, \code{\link{plot_bkwg}}.
  \item Estimation: \code{\link{fit_bkwg}}, \code{\link{bkwg_submodels}}.
  \item Inference: \code{\link{gof_bkwg}}, \code{\link{lrt_bkwg}},
    \code{\link{wald_bkwg}}, \code{\link{compare_bkwg}}.
}
}
\author{
Mohamad Yamen AL Mohamad \email{yamenmohamad@tabrizu.ac.ir}
}
\references{
Handique, L., Chakraborty, S. and Ali, M. M. (2017). Beta generated
Kumaraswamy-G family of distributions. \emph{Pakistan Journal of Statistics},
\bold{33}(6), 467--490. Preprint: \doi{10.48550/arXiv.1603.00634}.

Cordeiro, G. M. and de Castro, M. (2011). A new family of generalized
distributions. \emph{Journal of Statistical Computation and Simulation},
\bold{81}(7), 883--898.

Eugene, N., Lee, C. and Famoye, F. (2002). Beta-normal distribution and its
applications. \emph{Communications in Statistics -- Theory and Methods},
\bold{31}(4), 497--512.
}
\keyword{package}
\keyword{distribution}
