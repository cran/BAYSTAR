\name{TAR.thres}
\alias{TAR.thres}
\title{ To draw a threshold value. }
\description{
  The prior for the threshold parameter \eqn{thres}, follows a uniform prior on a range (l,u), where l and u can be set as relevant percentiles of the observed threshold variable. This prior could be considered to correspond to an empirical Bayes approach, rather than a fully Bayesian one.  The posterior distribution of \eqn{thres} is not of a standard distributional form, thus requiring us to use the Metropolis-Hastings (MH) method to achieve the desired sample for \eqn{thres}.

 }
\synopsis{
TAR.thres(ay, p1, p2, ph.1, ph.2, sig.1, sig.2, lagd, thres, step.r = 0.02, bound, lagp1, lagp2, constant = 1, thresVar)
}
\usage{
TAR.thres(ay, p1, p2, ph.1, ph.2, sig.1, sig.2, lagd, thres, step.r = 0.02, bound, lagp1, lagp2, constant = 1)

TAR.thres(ay, p1, p2, ph.1, ph.2, sig.1, sig.2, lagd, thres, step.r = 0.02, bound, lagp1, lagp2, constant = 1, thresVar)
}
\arguments{
  A list containing:
  \item{ay}{ The real data set. (input)    }
  \item{p1}{ Number of AR coefficients in regime one. }
  \item{p2}{ Number of AR coefficients in regime two. }
  \item{ph.1}{ The vector of AR parameters in regime one.}
  \item{ph.2}{ The vector of AR parameters in regime two.}
  \item{sig.1}{ The error terms of AR model in the regime one. }
  \item{sig.2}{ The error terms of AR model in the regime two. }
  \item{lagd}{ The delay lag parameter. }
  \item{thres}{ The threshold parameter. }
  \item{step.r}{ Step size of threshold variable for the MH algorithm are controlled the proposal variance. }
  \item{bound}{ The bound of threshold parameter. }
  \item{lagp1}{ The vector of non-zero autoregressive lags for the lower regime. (regime one); e.g. An AR model with p1=3, it could be non-zero lags 1,3, and 5 would set lagp1<-c(1,3,5).    }
  \item{lagp2}{ The vector of non-zero autoregressive lags for the upper regime. (regime two)   }
  \item{constant}{ Use the CONSTANT option to fit a model with/without a constant term (1/0). By default CONSTANT=1. }
  \item{thresVar}{ Exogenous threshold variable. (if missing, the series x is used) }
}

\author{ Cathy W.S. Chen, F.C. Liu }
\keyword{ methods }
