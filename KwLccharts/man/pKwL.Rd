\name{pKwL}
\alias{pKwL}
\title{Cumulative distribution function of the KwL distribution}
\description{CDF of the Kumaraswamy-Lindley distribution: \eqn{F(x)=1-[1-G(x)^a]^b}.}
\usage{pKwL(x, theta, a, b)}
\arguments{
  \item{x}{Non-negative numeric vector of quantiles.}
  \item{theta}{Lindley parameter, \eqn{\theta > 0}.}
  \item{a}{Kumaraswamy shape parameter, \eqn{a > 0}.}
  \item{b}{Kumaraswamy shape parameter, \eqn{b > 0}.}
}
\value{Numeric vector of cumulative probabilities.}
