\name{TAR.simu}
\alias{TAR.simu}
\title{ Simulated data from TAR model }
\description{
  To generate the simulated data from TAR(2;p1,p2) model.
}
\usage{
TAR.simu(nob, p1, p2, ph.1, ph.2, sig.1, sig.2, lagd, thres, lagp1, lagp2)
}
\arguments{
  \item{nob}{ Number of observations that we want to simulate. }
  \item{p1}{ Number of AR coefficient in regime one. }
  \item{p2}{ Number of AR coefficient in regime two. }
  \item{ph.1}{ The vector of AR parameters in regime one.}
  \item{ph.2}{ The vector of AR parameters in regime two.}
  \item{sig.1}{ The error terms in regime one. }
  \item{sig.2}{ The error terms in regime two. }  
  \item{lagd}{ The delay lag parameter. }
  \item{thres}{ The threshold parameter. }
  \item{lagp1}{ The vector of non-zero autoregressive lags for the lower regime. (regime one); e.g. An AR model with p1=3, it could be non-zero lags 1,3, and 5 would set lagp1<-c(1,3,5).    }
  \item{lagp2}{ The vector of non-zero autoregressive lags for the upper regime. (regime two)   }
}
\author{ Cathy W.S. Chen, Edward Lin }
\examples{
## Set the true values of all parameters
nob<- 2000                ## No. of observations
lagd<- 1                  ## delay lag of threshold variable
r<- 0.4                   ## r is the threshold value
sig.1<- 0.8; sig.2<- 0.5  ## variances of error distributions for two regimes
p1<- 2; p2<- 2            ## No. of covariate in two regimes
ph.1<- c(0.1,-0.4,0.3)    ## mean coefficients for regime 1
ph.2<- c(0.2,0.3,0.3)     ## mean coefficients for regime 2
lagp1<-1:2
lagp2<-1:2

yt<- TAR.simu(nob,p1,p2,ph.1,ph.2,sig.1,sig.2,lagd,r,lagp1,lagp2)
}
\keyword{ methods }
