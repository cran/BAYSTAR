\name{TAR.summary}
\alias{TAR.summary}
\title{ Calculate summary statistics for all parameters }
\description{
  A summary of the MCMC output can be obtained via the function \code{TAR.summary}.
  \code{TAR.summary} returns the posterior mean, median, standard deviation and
  the lower and upper bound of the 95\% Bayes posterior interval for all parameters,
  all botained from the sampling period only, after burn-in.
}
\usage{
TAR.summary(x, lagp1, lagp2, constant = 1)
}
\arguments{
  A list containing:
  \item{x}{ A matrix of the MCMC output of estimater parameters. }
  \item{lagp1}{ The vector of non-zero autoregressive lags for the lower regime. (regime one); e.g. An AR model with p1=3, it could be non-zero lags 1,3, and 5 would set lagp1<-c(1,3,5).    }
  \item{lagp2}{ The vector of non-zero autoregressive lags for the upper regime. (regime two)   }
  \item{constant}{ Use the CONSTANT option to fit a model with/without a constant term (1/0). By default CONSTANT=1. }
}
\author{ Cathy W.S. Chen, F.C. Liu }
\keyword{ methods }