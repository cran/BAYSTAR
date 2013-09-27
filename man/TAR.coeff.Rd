\name{TAR.coeff}
\alias{TAR.coeff}
\title{ Estimate AR coefficients }
\description{
  We assume a normal prior for the AR coefficients and draw AR coefficients from a multivariate normal posterior distribution. Parsimonious subset AR could be assigned in each regime in the BAYSTAR function rather than a full AR model.
}
\usage{
TAR.coeff(reg, ay, p1, p2, sig, lagd, 
          thres, mu0, v0, lagp1, lagp2, constant = 1, thresVar)
}
\arguments{
  A list containing:
  \item{reg}{ The regime is assigned. (equal to one or two)  }
  \item{ay}{ The real data set. (input)   }
  \item{p1}{ Number of AR coefficients in regime one. }
  \item{p2}{ Number of AR coefficients in regime two. }
  \item{sig}{ The error terms of TAR model. }
  \item{lagd}{ The delay lag parameter. }
  \item{thres}{ The threshold parameter. }
  \item{mu0}{ Mean vector of conditional prior distribution in mean equation. }
  \item{v0}{ Covariance matrix of conditional prior distribution in mean equation. }
  \item{lagp1}{ The vector of non-zero autoregressive lags for the lower regime. (regime one); e.g. An AR model with p1=3, it could be non-zero lags 1,3, and 5 would set lagp1<-c(1,3,5).    }
  \item{lagp2}{ The vector of non-zero autoregressive lags for the upper regime. (regime two)   }
  \item{constant}{ Use the CONSTANT option to fit a model with/without a constant term (1/0). By default CONSTANT=1. }
  \item{thresVar}{ Exogenous threshold variable. (if missing, the self series are used)}
}

\author{ Cathy W.S. Chen, F.C. Liu }

\keyword{ methods }
