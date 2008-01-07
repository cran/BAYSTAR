`BAYSTAR` <-
function(x,lagp1,lagp2,nIteration,nBurnin,constant=1,differ=0,d0=3,step.thv,thresVar) {
##Time.initial<-Sys.time()
## Initialize

p1<- length(lagp1); p2<- length(lagp2)            ## No. of covariate in two regimes
nx<-length(x)
if (differ ==1){
yt<-x[2:nx]-x[2:nx-1]   }
else yt <- x

nob<- length(yt)

if (!missing(thresVar)){
    if (length(thresVar) > nob ){
        zt <- thresVar[1:nob]}
    else zt<-thresVar
}
else zt<-yt         
## Set initial values
phi.1<- rep(0.05,p1)
phi.2<- rep(0.05,p2)
sigma.1<- 0.2
sigma.2<- 0.2
lagd<- 3
thres<- median(zt)
accept.r<- 0
sum.r<- 0

## MSE of fitting an AR(p1) model
ar.mse<- ar(yt,aic=FALSE, order.max=p1)

## Sets for the hyper-parameters
mu01<- matrix(0,nrow=p1+constant,ncol=1)
v01<- diag(0.1,p1+constant)
mu02<- matrix(0,nrow=p2+constant,ncol=1)
v02<- diag(0.1,p2+constant)
v<- 3
lambda<- ar.mse$var.pred/3
bound.thv<- c(quantile(zt,0.25),quantile(zt,0.75))





## Initialize a matrix for saving all iterative estimates
if(constant==1){
par.set<- matrix(NA,nrow=nIteration,ncol=(length(c(phi.1,phi.2,sigma.1,sigma.2,lagd,thres))+2))}
else{
par.set<- matrix(NA,nrow=nIteration,ncol=length(c(phi.1,phi.2,sigma.1,sigma.2,lagd,thres)))}
## Start of MCMC sampling
for (igb in 1:nIteration){
if (!missing(thresVar)){
phi.1<- TAR.coeff(1,yt,p1,p2,sigma.1,lagd,thres,mu01,v01,lagp1,lagp2,constant=constant,zt)## Draw phi.1 from a multivariate normal distribution
phi.2<- TAR.coeff(2,yt,p1,p2,sigma.2,lagd,thres,mu02,v02,lagp1,lagp2,constant=constant,zt)## Draw phi.2 from a multivariate normal distribution
sigma.1<- TAR.sigma(1,yt,thres,lagd,p1,p2,phi.1,v,lambda,lagp1,lagp2,constant=constant,zt)      ## Draw sigma.1 from an Inverse-Gamma distribution  ## v and lambda are the hyper-parameters of the Gamma prior
sigma.2<- TAR.sigma(2,yt,thres,lagd,p1,p2,phi.2,v,lambda,lagp1,lagp2,constant=constant,zt)      ## Draw sigma.2 from a Inverse-Gamma distribution
lagd<- TAR.lagd(yt,p1,p2,phi.1,phi.2,sigma.1,sigma.2,thres,lagp1,lagp2,constant=constant,d0,zt)    ## Draw lagd from a multinomial distribution
thresholdt<- TAR.thres(yt,p1,p2,phi.1,phi.2,sigma.1,sigma.2,lagd,thres,step.r=step.thv,bound.thv,lagp1,lagp2,constant=constant,zt)            ## Draw thresholdt by the MH algorithm
}
else{
phi.1<- TAR.coeff(1,yt,p1,p2,sigma.1,lagd,thres,mu01,v01,lagp1,lagp2,constant=constant)## Draw phi.1 from a multivariate normal distribution
phi.2<- TAR.coeff(2,yt,p1,p2,sigma.2,lagd,thres,mu02,v02,lagp1,lagp2,constant=constant)## Draw phi.2 from a multivariate normal distribution
sigma.1<- TAR.sigma(1,yt,thres,lagd,p1,p2,phi.1,v,lambda,lagp1,lagp2,constant=constant)      ## Draw sigma.1 from an Inverse-Gamma distribution  ## v and lambda are the hyper-parameters of the Gamma prior
sigma.2<- TAR.sigma(2,yt,thres,lagd,p1,p2,phi.2,v,lambda,lagp1,lagp2,constant=constant)      ## Draw sigma.2 from a Inverse-Gamma distribution
lagd<- TAR.lagd(yt,p1,p2,phi.1,phi.2,sigma.1,sigma.2,thres,lagp1,lagp2,constant=constant,d0)    ## Draw lagd from a multinomial distribution
thresholdt<- TAR.thres(yt,p1,p2,phi.1,phi.2,sigma.1,sigma.2,lagd,thres,step.r=step.thv,bound.thv,lagp1,lagp2,constant=constant)            ## Draw thresholdt by the MH algorithm
}
sum.r<- sum.r+thresholdt[1]   ## Count the number of acceptance
thres<- thresholdt[2]         ## Save i-th iterated threshold value

if(constant==1){
c.mean<- c(phi.1[1]/(1-sum(phi.1)+phi.1[1]),phi.2[1]/(1-sum(phi.2)+phi.2[1]))
par.set[igb,]<-c(phi.1,phi.2,sigma.1,sigma.2,thres,c.mean,lagd)
} ## Compute the unconditional means for each regime
else {par.set[igb,]<-c(phi.1,phi.2,sigma.1,sigma.2,thres,lagd)
}
      ## Save all iterated estimates of parameters
ncol0<-ncol(par.set)

## Print out for monitoring the estimations of every 1000 iterate
if(igb%%1000==0){
cat("iteration = ",igb,"\n")
cat("regime 1 = ",round(phi.1,4),"\n")
cat("regime 2 = ",round(phi.2,4),"\n")
cat("sigma 1  = ",round(sigma.1,4),"\n")
cat("sigma 2  = ",round(sigma.2,4),"\n")
cat("r        = ",round(thres,4),"\n")
accept.r<- (sum.r/igb)*100
cat("acceptance rate of r = ", round(accept.r,4),"%", "\n")

## Make a frequency table of delay lag
lag.freq<- rep(0,d0)
for(i in 1:d0){
lag.freq[i]<- sum(par.set[1:igb,ncol0]==i)
}
#lag.freq[1:length(table(par.set[,ncol0]))]<- table(par.set[,ncol0]) ## Frequency table of delay lag
lag.freq<- t(matrix(lag.freq,dimnames=list(c(as.character(1:d0)),c("Freq"))))
cat("Lag choice : ", "\n")
print(lag.freq)
cat("------------","\n")
}
} ## End of MCMC sampling
## Summarize the collected MCMC estimates
mcmc.stat<- TAR.summary(par.set[(nBurnin+1):nIteration,1:(ncol0-1)],lagp1,lagp2,constant=constant)
print(round(mcmc.stat,4))
## Calculate the highest posterior probability of delay lag
lag.y<- c(1:d0)
lag.d<- lag.y[lag.freq==max(lag.freq)]
cat("Lag choice : ", "\n")
print(lag.freq)
cat("------------","\n")
cat("The highest posterior prob. of lag at : ",lag.d,"\n")

## Calculate the residual for TAR model
maxd<-max(lagp1,lagp2)
if (constant == 1){
con.1<-mcmc.stat[1,1]
par.1<-mcmc.stat[2:(p1+1),1]
con.2<-mcmc.stat[p1+2,1]
par.2<-mcmc.stat[(p1+2+1):(p1+p2+2),1]
thv  <-mcmc.stat[p1+p2+3,1]
}else{par.1<-mcmc.stat[1:p1,1]
par.2<-mcmc.stat[(p1+1):(p1+p2),1]
thv  <-mcmc.stat[p1+p2+1,1]
}
residual<-rep(NA,nob-maxd)
for (t in (maxd+1):nob){
if (constant == 1){
 if ( yt[t-lag.d] <= thv){
  residual[t-maxd]<- yt[t] - sum(con.1,(par.1 * yt[t-lagp1]))
 }
 else{
  residual[t-maxd]<- yt[t] - sum(con.2,(par.2 * yt[t-lagp2]))
 }
}
else{
 if ( yt[t-lag.d] <= thv){
  residual[t-maxd]<- yt[t] - sum(par.1 * yt[t-lagp1])
 }
 else{
  residual[t-maxd]<- yt[t] - sum(par.2 * yt[t-lagp2])
 }
}
}
tar<-list(mcmc=par.set,coef=round(mcmc.stat,4),residual=residual,lagd=lag.d)
return(tar)
##Sys.time()-Time.initial
}

