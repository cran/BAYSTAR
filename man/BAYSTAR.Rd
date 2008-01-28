\name{BAYSTAR}
\alias{BAYSTAR}
\title{ Threshold Autoregressive model: Bayesian approach }
\description{
  Bayesian estimation and inference for two-regime TAR model, as well as monitoring MCMC convergence. One may want to allow for higher-order AR models in the different regimes. Parsimonious subset AR could be assigned in each regime in the BAYSTAR function rather than a full AR model (i.e. Autoregressive order could be not a continuous series).
}
\synopsis{
BAYSTAR(x, lagp1, lagp2, nIteration, nBurnin, constant = 1, differ = 0, d0 = 3, step.thv, thresVar)
}
\usage{
BAYSTAR(x, lagp1, lagp2, nIteration, nBurnin, constant = 1, differ = 0, d0 = 3, step.thv)

BAYSTAR(x, lagp1, lagp2, nIteration, nBurnin, constant = 1, differ = 0, d0 = 3, step.thv, thresVar)
}

\arguments{
  \item{x}{ Time series. }
  \item{lagp1}{ The vector of non-zero autoregressive lags for the lower regime. (regime one); e.g. An AR model with p1=3, it could be non-zero lags 1,3, and 5 would set lagp1<-c(1,3,5).    }
  \item{lagp2}{ The vector of non-zero autoregressive lags for the upper regime. (regime two)   }
  \item{nIteration}{ Total MCMC iterations. }
  \item{nBurnin}{ Burn-in iterations. }
  \item{constant}{ Use the CONSTANT option to fit a model with/without a constant term (1/0). By default CONSTANT=1. }
  \item{differ}{ Take the first difference. (default = 0)  }
  \item{d0}{The option of a set maximum delay. (default = 3)}
  \item{step.thv}{ Step size of threshold variable for the MH algorithm are controlled the proposal variance. }
  \item{thresVar}{ Exogenous threshold variable. (if missing, the series x is used.) }
}
\details{
  Given the maximum AR orders p1 and p2, the two-regime SETAR(2:p1;p2) model is specified as:
\deqn{x_{t} = ( \phi _0^{(1)} + \phi _1^{(1)} x_{t - 1} + \ldots + \phi _{p_1 }^{(1)} x_{t - p_1 } + a_t^{(1)} ) I( z_{t-d} <= th) + ( \phi _0^{(2)}  + \phi _1^{(2)} x_{t - 1} + \ldots + \phi _{p_2 }^{(2)} x_{t - p_2 } + a_t^{(2)} I( z_{t-d} > th) }

where \eqn{th} is the threshold parameter driving the regime-switching behavior; \eqn{z_{t}} is the threshold variable; \eqn{d} is the threshold lag of the model; and the error term \eqn{a_t^{(j)}} in regime \eqn{j, (j=1,2)} is assumed to be an i.i.d. Gaussian white noise process with mean zero and variance \eqn{sigma_j^2, j=1,2}. I[A] is an indicator function with I[A]=1 if the event A occurs and I[A]=0 otherwise. One may want to allow parsimonious subset AR model in each regime rather than a full AR model.
}

\value{
A list with output containing the following components:
  \item{mcmc }{the mcmc results of all parameters (including burn-in).}
  \item{coef }{summarize the collected MCMC estimates after burn-in. Including all parameters of TAR model.}
  \item{residual }{the residuals from the fitting model by \code{coef}.}
  \item{lagd }{The mode of time delay of threshold variable in MCMC iteration.}
}

\author{ Cathy W. S. Chen, Edward M.H. Lin, F.C. Liu, and Richard Gerlach  }
\examples{
data(unemployrate)
x<- unemployrate
lagp1<-c(2,3,4,10,12)
lagp2<-c(2,3,12)
## Total MCMC iterations and burn-in iterations
nIterations<- 10000
nBurnin<- 2000
## Step size for the MH algorithm
step.thv<- 2.0
y<-BAYSTAR(x,lagp1,lagp2,nIterations,nBurnin,constant=0,differ=1,d0=3,step.thv=step.thv)
}

\keyword{ methods }
