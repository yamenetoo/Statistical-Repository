\name{ross_ellipse}
\alias{ross_ellipse}
\docType{data}
\title{Ellipse Shape Data (Table 1.2)}
\description{
Twelve values of the ratio \eqn{x = a/b} of the axes of an ellipse,
classified into two classes, from the inductive-reasoning example of
Ross (2000, Table 1.2).
}
\usage{data(ross_ellipse)}
\format{
A data frame with 12 rows and columns \code{x} (the ratio \eqn{a/b}) and
\code{class} (1 or 2).
}
\source{
Ross, T. J. (2000). Membership functions, fuzzification and defuzzification.
In \emph{Fuzzy Systems in Medicine}, 49--78. Springer-Verlag.
}
\examples{
data(ross_ellipse)
entropy_partition(ross_ellipse$x, ross_ellipse$class, levels = 2)
}
\keyword{datasets}
