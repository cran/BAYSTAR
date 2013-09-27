\name{TAR.lagd}
\alias{TAR.lagd}
\title{ Identification of lag order of threshold variable }
\description{
  The delay d has a discrete uniform prior over the integers: 1,2,..., d0, where d0 is a set maximum delay.  We draw the delay lag of threshold variable from a multinomial distribution.
}
\usage{
TAR.lagd(ay, p1, p2, ph.1, ph.2, sig.1, sig.2, 
         thres, lagp1, lagp2, constant = 1, d0, thresVar)
}
\arguments{
  A list containing:
  \item{ay}{ The real data set. (input) }
  \item{p1}{ Number of AR coefficients in regime one.  }
  \item{p2}{ Number of AR coefficients in regime two.  }
  \item{ph.1}{ The vector of AR parameters in regime one.}
  \item{ph.2}{ The vector of AR parameters in regime two.}
  \item{sig.1}{ The error terms of AR model in the regime one. }
  \item{sig.2}{ The error terms of AR model in the regime two. }  
  \item{thres}{ The threshold parameter. }
  \item{lagp1}{ The vector of non-zero autoregressive lags for the lower regime. (regime one); e.g. An AR model with p1=3, it could be non-zero lags 1,3, and 5 would set lagp1<-c(1,3,5).    }
  \item{lagp2}{ The vector of non-zero autoregressive lags for the upper regime. (regime two)   }
  \item{constant}{ Use the CONSTANT option to fit a model with/without a constant term (1/0). By default CONSTANT=1. }
  \item{d0}{The maximum delay lag could be selected.} 
  \item{thresVar}{ Exogenous threshold variable. (if missing, the series x is used) }
}
\author{ Cathy W.S. Chen, Edward Lin }

\keyword{ methods }

