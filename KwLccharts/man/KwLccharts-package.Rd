\name{KwLccharts-package}
\alias{KwLccharts-package}
\alias{KwLccharts}
\docType{package}
\title{Bootstrap Percentile Control Charts for the Kumaraswamy-Lindley Distribution}
\description{
Provides the density, cumulative distribution, closed-form quantile (via
the Lambert W function) and random-generation functions for the
Kumaraswamy-Lindley (KwL) distribution, together with maximum-likelihood
fitting (\code{\link{fit_KwL}}) and a parametric bootstrap percentile
control chart (\code{\link{KwLBootCC}}, \code{\link{KwL_ARL}}), following
the methodology of Al-Lami, Jabbari Khamnei and Heydari (2025, PLOS ONE)
for the Lindley-Geometric distribution, adapted here to the three-parameter
Kumaraswamy-Lindley family of Cordeiro and de Castro (2011).
}
\details{
Main functions: \code{\link{dKwL}}, \code{\link{pKwL}}, \code{\link{qKwL}},
\code{\link{rKwL}}, \code{\link{gLindley}}, \code{\link{GLindley}},
\code{\link{QLindley}}, \code{\link{fit_KwL}}, \code{\link{KwLBootCC}},
\code{\link{KwL_ARL}}.
}
\author{
Mohamad Yamen Al-Mohamad
}
