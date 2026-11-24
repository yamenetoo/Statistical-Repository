\name{ross_colors}
\alias{ross_colors}
\docType{data}
\title{Pairwise Colour Preferences (Table 1.1)}
\description{
Pairwise preferences of 1000 people among five colours, from the
rank-ordering example of Ross (2000, Table 1.1).
}
\usage{data(ross_colors)}
\format{
A 5 x 5 numeric matrix with row and column names red, orange, yellow, green
and blue: entry \code{[i, j]} is the number of people who preferred colour
\code{i} to colour \code{j}; the diagonal is \code{NA}.
}
\source{
Ross, T. J. (2000). Membership functions, fuzzification and defuzzification.
In \emph{Fuzzy Systems in Medicine}, 49--78. Springer-Verlag.
}
\examples{
data(ross_colors)
rank_ordering(ross_colors)
}
\keyword{datasets}
