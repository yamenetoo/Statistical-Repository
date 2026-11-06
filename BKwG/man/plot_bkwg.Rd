\name{plot_bkwg}
\alias{plot_bkwg}
\title{Plot the Density, CDF, Survival or Hazard of a BKw-G Distribution}
\description{
Draws the pdf, cdf, survival or hazard function of a BKw-G distribution,
useful for exploring the shapes of the family.
}
\usage{
plot_bkwg(a = 1, b = 1, alpha = 1, beta = 1, baseline = "weibull",
          par = NULL, what = c("pdf", "cdf", "survival", "hazard"),
          from = NULL, to = NULL, add = FALSE, n = 401, ...)
}
\arguments{
  \item{a, b, alpha, beta}{positive shape parameters.}
  \item{baseline, par}{baseline distribution and its parameters, as in
    \code{\link{dbkwg}}.}
  \item{what}{the function to draw.}
  \item{from, to}{plotting range; by default the 0.001 and 0.999 quantiles.}
  \item{add}{logical; add to the current plot.}
  \item{n}{number of evaluation points.}
  \item{\dots}{graphical parameters passed to \code{plot} or \code{lines}.}
}
\value{Invisibly, a data frame with the plotted coordinates.}
\seealso{\code{\link{dbkwg}}}
\examples{
plot_bkwg(a = 2, b = 1, alpha = 1.5, beta = 2, par = c(2, 1),
          to = 3, ylim = c(0, 1.6))
plot_bkwg(a = 0.5, b = 3, alpha = 2, beta = 0.5, par = c(2, 1),
          from = 0, to = 3, add = TRUE, col = "red")
plot_bkwg(a = 2, b = 0.4, alpha = 0.5, beta = 3, par = c(1, 1),
          what = "hazard", to = 4)
}
\keyword{hplot}
