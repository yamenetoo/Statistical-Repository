\name{plot_fuzzy_sets}
\alias{plot_fuzzy_sets}
\alias{plot.fuzzyset}
\alias{lines.fuzzyset}
\title{Plot Fuzzy Sets}
\description{Plot one or several membership functions.}
\usage{
plot_fuzzy_sets(sets, col = NULL, lwd = 2, lty = 1, xlim = NULL,
                ylim = c(0, 1.05), xlab = "x", ylab = expression(mu(x)),
                main = "", legend = TRUE, legend_pos = "topright", \dots)

\method{plot}{fuzzyset}(x, xlim = NULL, ylim = c(0, 1.05), xlab = "x",
     ylab = expression(mu(x)), main = if (nzchar(x$name)) x$name else "",
     col = "black", lwd = 2, \dots)
\method{lines}{fuzzyset}(x, col = "black", lwd = 2, \dots)
}
\arguments{
  \item{sets}{a list of \code{"fuzzyset"} objects.}
  \item{x}{a \code{"fuzzyset"} object.}
  \item{col, lwd, lty}{line colours, widths and types.}
  \item{xlim, ylim, xlab, ylab, main}{usual plot arguments.}
  \item{legend}{draw a legend with the set names?}
  \item{legend_pos}{position of the legend.}
  \item{\dots}{further arguments passed to \code{plot} or \code{lines}.}
}
\value{The input, invisibly.}
\examples{
plot(fs_trapezoidal(1, 3, 5, 8, name = "trapezoid"))
plot_fuzzy_sets(ross2000_sets("height"))      # Figure 1.1
}
\keyword{hplot}
