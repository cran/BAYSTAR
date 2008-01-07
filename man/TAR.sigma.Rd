\name{TAR.sigma}
\alias{TAR.sigma}
\title{ To draw the variance of error distribution. }
\description{
  We employ a conjugate prior, Inverse-Gamma distribution, for sigma squared in regime j, j=1,2. To draw the variance of error distribution from an Inverse-Gamma posterior distribution.
}
\synopsis{
TAR.sigma(reg, ay, thres, lagd, p1, p2, ph, v, lambda, lagp1, lagp2, constant = 1, thresVar)
}
\usage{
TAR.sigma(reg, ay, thres, lagd, p1, p2, ph, v, lambda, lagp1, lagp2, constant = 1)

TAR.sigma(reg, thres, lagd, p1, p2, ph, ay, v, lambda, lagp1, lagp2, constant = 1, thresVar)
}

\arguments{
  A list containing:
  \item{reg}{ The regime is assigned. (equal to one or two)  }
  \item{thres}{ The threshold parameter. }
  \item{lagd}{ The delay lag parameter. }
  \item{p1}{ Number of AR coefficient in regime one. }
  \item{p2}{ Number of AR coefficient in regime two. }
  \item{ph}{ The vector of AR parameters in regime \code{reg}.}
  \item{ay}{ The real data set. (input)   }
  \item{v, lambda}{ The hyper-parameter of Inverse Gamma distribution for priors of variance. (i.e. IG(v/2,lambda/2)) }
  \item{lagp1}{ The vector of non-zero autoregressive lags for the lower regime. (regime one); e.g. An AR model with p1=3, it could be non-zero lags 1,3, and 5 would set lagp1<-c(1,3,5).    }
  \item{lagp2}{ The vector of non-zero autoregressive lags for the upper regime. (regime two)   }
  \item{constant}{ Use the CONSTANT option to fit a model with/without a constant term (1/0). By default CONSTANT=1. }
  \item{thresVar}{ Exogenous threshold variable. (if missing, the series x is used) }
}

\author{ Cathy W.S. Chen, Edward Lin }

\keyword{ methods }

