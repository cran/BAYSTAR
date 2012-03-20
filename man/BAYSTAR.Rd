\name{BAYSTAR}
\alias{BAYSTAR}
\title{ Threshold Autoregressive model: Bayesian approach }
\description{
Bayesian estimation and one-step-ahead forecasting for two-regime
TAR model, as well as monitoring MCMC convergence. One may want to
allow for higher-order AR models in the different regimes.
Parsimonious subset AR could be assigned in each regime in the
BAYSTAR function rather than a full AR model (i.e. the
autoregressive orders could be not a continuous series).}
\synopsis{
BAYSTAR(x, lagp1, lagp2, Iteration, Burnin, constant, d0, step.thv, thresVar, mu01, v01, mu02, v02, v0, lambda0, refresh)
}
\usage{

BAYSTAR(x, lagp1, lagp2, Iteration, Burnin, constant, d0, step.thv, thresVar, mu01, v01, mu02, v02, v0, lambda0, refresh)

}

\arguments{
  \item{x}{ A vector of time series. }
  
  \item{lagp1}{ A vector of non-zero autoregressive lags for the lower regime (regime one).
                For example, an AR model with p1=3 in lags 1,3, and 5 would be set
                as \code{lagp1<-c(1,3,5)}.}
  
  \item{lagp2}{ A vector of non-zero autoregressive lags for the upper regime (regime two).}
  
  \item{Iteration}{ The number of MCMC iterations.}
  
  \item{Burnin}{ The number of burn-in iterations for the sampler.}
  
  \item{constant}{ The intercepts include in the model for each regime, if \code{constant}=1.
                   Otherwise, if \code{constant}=0. (Default: \code{constant}=1)}
  
  \item{d0}{ The maximum delay lag considered. (Default: \code{d0} = 3)}
  
  \item{step.thv}{ Step size of tuning parameter for the Metropolis-Hasting algorithm. }
  
  \item{thresVar}{ A vector of time series for the threshold variable. (if missing, the series x is used.) }
  
  \item{mu01}{ The prior mean of \eqn{\phi}{phi} in regime one.
               This setting can be a scalar or a column vector with dimension equal
               to the number of \eqn{\phi}{phi}. If this sets a scalar value, then
               the prior mean for all of \eqn{\phi}{phi} are this value. (Default:
               a vector of zeros)}
  
  \item{v01}{ The prior covariance matrix of \eqn{\phi}{phi} in regime one.
              This setting can either be a scalar or a square matrix with
              dimensions equal to the number of \eqn{\phi}{phi}. If this sets a
              scalar value, then prior covariance matrix of \eqn{\phi}{phi} is
              that value times an identity matrix.  (Default: a diagonal matrix
              are set to 0.1)}
  
  \item{mu02}{ The prior mean of \eqn{\phi}{phi} in regime two.
               This setting can be a scalar or a column vector with dimension equal
               to the number of \eqn{\phi}{phi}. If this sets a scalar value, then
               the prior mean for all of \eqn{\phi}{phi} are this value. (Default:
               a vector of zeros)}
  
  \item{v02}{ The prior covariance matrix of \eqn{\phi}{phi} in regime two.
              This setting can either be a scalar or a square matrix with
              dimensions equal to the number of \eqn{\phi}{phi}. If this sets a
              scalar value, then prior covariance matrix of \eqn{\phi}{phi} is
              that value times an identity matrix.  (Default: a diagonal matrix
              are set to 0.1)}
  
  \item{v0}{ \code{v0}/2 is the shape parameter for Inverse-Gamma prior of \eqn{\sigma^2}{sigma^2}.
             (Default: \code{v0} = 3)}
  
  \item{lambda0}{ \code{lambda0}*\code{v0}/2 is the scale parameter for Inverse-Gamma prior of
                  \eqn{\sigma^2}{sigma^2}. (Default: \code{lambda0} = the residual
                  mean squared error of fitting an AR(p1) model to the data.)}
  
  \item{refresh}{ Each \code{refresh} iteration for monitoring MCMC output. (Default: \code{refresh}=1000)}
  
 
}
\details{
  Given the maximum AR orders p1 and p2, the two-regime SETAR(2:p1;p2) model is specified as:
\deqn{x_{t} = ( \phi _0^{(1)} + \phi _1^{(1)} x_{t - 1} + \ldots + \phi _{p1 }^{(1)} x_{t - p1 } + 
      a_t^{(1)} ) I( z_{t-d} <= th) + ( \phi _0^{(2)}  + \phi _1^{(2)} x_{t - 1} + \ldots + 
      \phi _{p2 }^{(2)} x_{t - p2 } + a_t^{(2)} I( z_{t-d} > th) }

where \eqn{th}{} is the threshold value for regime switching;
\eqn{z_{t}}{} is the threshold variable; \eqn{d}{} is the delay lag
of threshold variable; and the error term \eqn{a_t^{(j)}}{}, \eqn{j,
(j=1,2)}{}, for each regime is assumed to be an i.i.d. Gaussian
white noise process with mean zero and variance \eqn{sigma_j^2,
j=1,2}{}. I(A) is an indicator function. Event A will occur if
I(A)=1 and otherwise if I(A)=0. One may want to allow parsimonious
subset AR model in each regime rather than a full AR model.
}

\value{
A list of output with containing the following components:
  \item{mcmc }{All MCMC iterations.}
  \item{posterior }{The initial \code{Burnin} iterations are discarded as a burn-in sample, the final sample of (\code{Iteration-Burnin}) iterates is used for posterior inference.}
  \item{coef }{Summary Statistics of parameter estimation based on the final sample of (\code{Iteration-Burnin}) iterates.}
  \item{residual }{Residuals from the estimated model.}
  \item{lagd }{The mode of time delay lag of the threshold variable.}
  \item{DIC }{The deviance information criterion (DIC); a Bayesian method for model comparison (Spiegelhalter et al, 2002)}
}

\author{ Cathy W. S. Chen, Edward M.H. Lin, F.C. Liu, and Richard Gerlach  }
\examples{
library(BAYSTAR)
data(unemployrate)
x  <- unemployrate
nx <- length(x)
dif.x <- x[2:nx]-x[2:nx-1]
lagp1 <- c(2,3,4,10,12)
lagp2 <- c(2,3,12)

## Total MCMC iterations and burn-in iterations
Iteration <- 10000
Burnin    <- 2000

## A RW (random walk) MH algorithm is used in simulating the threshold value
## Step size for the RW MH
step.thv<- 2.5

out <- BAYSTAR(dif.x,lagp1,lagp2,Iteration,Burnin,constant=0,step.thv=step.thv)

library(coda)
geweke.diag(out$posterior, frac1=0.1, frac2=0.5)
}

\keyword{ methods }
