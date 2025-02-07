### 
###  Simulate 1-D SPP 
### 

S.domain=c(0,1)

n.win=3
tot.area=S.domain[2]-S.domain[1]
tot.win.area=.5*tot.area
win.area=tot.win.area/n.win
tot.nonwin.area=tot.area-tot.win.area
nonwin.area=tot.nonwin.area/(n.win+1)

S.subdomain.mat=matrix(cumsum(rep(c(nonwin.area,win.area),n.win)),,2,byrow=TRUE)

m=10000 # number of quadrature points
s.full=seq(S.domain[1],S.domain[2],,m)
X.full=matrix(0,m,2)
X.full[,1]=s.full
X.full[,2]=sin(20*s.full)
ds=X.full[2,1]-X.full[1,1]

beta=c(2,1)
beta.0=4
lam.full=exp(beta.0+X.full%*%beta)
lam.max=max(lam.full)

M=rpois(1,lam.max) 
s.superpop=runif(M,S.domain[1],S.domain[2]) # potential starting points
X.superpop=cbind(s.superpop,sin(20*s.superpop))
lam.superpop=exp(beta.0+X.superpop%*%beta)

tot.idx=rbinom(M,1,lam.superpop/lam.max)==1
s.tot=s.superpop[tot.idx] # total observed points 
X.tot=X.superpop[tot.idx,] 
N=length(s.tot)

###
###  Get windowed data 
###

win.idx=NULL
win.full.idx=NULL
for(j in 1:n.win){
	tmp.win.idx=(1:N)[(s.tot>S.subdomain.mat[j,1])&(s.tot<S.subdomain.mat[j,2])]
  win.idx=c(win.idx,tmp.win.idx)
	tmp.full.idx=(1:m)[(s.full>S.subdomain.mat[j,1])&(s.full<S.subdomain.mat[j,2])]
  win.full.idx=c(win.full.idx,tmp.full.idx)
}
n=length(win.idx)

s.win=s.tot[win.idx]
X.win=X.tot[win.idx,]

X.win.full=X.full[win.full.idx,]
X.nowin.full=X.full[-win.full.idx,]
win.area=sum(S.subdomain.mat[,2]-S.subdomain.mat[,1])

###
###  Make Simulation Plots
###

layout(matrix(1:3,3,1))
plot(function(x){exp(beta.0+beta[1]*x+beta[2]*sin(20*x))},xlab="location",ylab=bquote(lambda),from=S.domain[1],to=S.domain[2],n=1000)
apply(S.subdomain.mat,1,function(x){rect(x[1],par("usr")[3],x[2],par("usr")[4],col=rgb(0,0,0,.05),border=NA)})
rug(s.superpop,.05,col=rgb(0,0,0,.2));rug(s.tot,.05,lwd=1.5,col=rgb(1,0,0,.5));rug(s.win,.05,lwd=1.5,col=3)
plot(function(x){x},xlab="location",ylab=bquote(x[1]),from=S.domain[1],to=S.domain[2],n=1000)
apply(S.subdomain.mat,1,function(x){rect(x[1],par("usr")[3],x[2],par("usr")[4],col=rgb(0,0,0,.05),border=NA)})
rug(s.tot,.05,lwd=1.5,col=rgb(1,0,0,.5));rug(s.win,.05,lwd=1.5,col=3)
plot(function(x){sin(20*x)},xlab="location",ylab=bquote(x[2]),from=S.domain[1],to=S.domain[2],n=1000)
apply(S.subdomain.mat,1,function(x){rect(x[1],par("usr")[3],x[2],par("usr")[4],col=rgb(0,0,0,.05),border=NA)})
rug(s.tot,.05,lwd=1.5,col=rgb(1,0,0,.5));rug(s.win,.05,lwd=1.5,col=3)

###
###  Fit SPP w/ comp. lik. to full data
###

n.mcmc=100000
source("spp.comp.mcmc.R")
out.comp.full=spp.comp.mcmc(s.win,X.win,X.win.full,ds,win.area,n.mcmc)

layout(matrix(1:2,2,1))
plot(out.comp.full$beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(out.comp.full$beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

###
### Fit SPP w/ comp. lik. w/ ESN
###
source(here("GlacierBay_Code", "spp.comp.ESN.mcmc.R"))
out.comp.esn=spp.comp.ESN.mcmc(s.win,cbind(rep(1, nrow(X.win)), X.win),
                               cbind(rep(1, nrow(X.win.full)), X.win.full),
                               5,ds,n.mcmc,0.1)

###
###  Fit SPP w/ cond. lik. to full data
###

n.mcmc=100000
source("spp.cond.mcmc.R")
out.cond.full=spp.cond.mcmc(s.win,X.win,X.win.full,ds,n.mcmc)

layout(matrix(1:2,2,1))
plot(out.cond.full$beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(out.cond.full$beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

###
###  Fit SPP using cond. output using 2nd stage MCMC 
###

source("spp.stg2.mcmc.R")
out.cond.2.full=spp.stg2.mcmc(out.cond.full)

layout(matrix(1:2,2,1))
plot(out.cond.2.full$beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(out.cond.2.full$beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

###
###  Compare marginal posteriors 
###

layout(matrix(1:3,1,3))
hist(out.comp.full$beta.0.save,prob=TRUE,breaks=60,main="",xlab=bquote(beta[0]))
lines(density(out.cond.full$beta.0.save,n=1000),col=2,lwd=2)
lines(density(out.cond.2.full$beta.0.save,n=1000),col=3,lwd=2)
hist(out.comp.full$beta.save[1,],prob=TRUE,breaks=60,main="",xlab=bquote(beta[1]))
lines(density(out.cond.full$beta.save[1,],n=1000),col=2,lwd=2)
lines(density(out.cond.2.full$beta.save[1,],n=1000),col=3,lwd=2)
hist(out.comp.full$beta.save[2,],prob=TRUE,breaks=60,main="",xlab=bquote(beta[2]))
lines(density(out.cond.full$beta.save[2,],n=1000),col=2,lwd=2)
lines(density(out.cond.2.full$beta.save[2,],n=1000),col=3,lwd=2)

###
###  Posterior for N 
###

N.save=rep(0,n.mcmc)

for(k in 1:n.mcmc){
	if(k%%10000==0){cat(k," ")}
	beta.0.tmp=out.cond.2.full$beta.0.save[k]
	beta.tmp=out.cond.2.full$beta.save[,k]
  lam.nowin.int=sum(exp(log(ds)+beta.0.tmp+X.nowin.full%*%beta.tmp))
  N.save[k]=n+rpois(1,lam.nowin.int)
};cat("\n")

plot(N.save,type="l")
abline(h=N,col=rgb(0,1,0,.8),lty=2,lwd=2)

hist(N.save,breaks=50,prob=TRUE,main="",xlab="N")
abline(v=N,col=rgb(0,1,0,.8),lty=2,lwd=2)

###
###  PPD of lambda for Full Area 
###

idx.sm=seq(1,m,10)
m.sm=length(idx.sm)
s.sm=s.full[idx.sm]
X.sm=X.full[idx.sm,]
lam.save=matrix(0,m.sm,n.mcmc)
for(k in 1:n.mcmc){
	if(k%%10000==0){cat(k," ")}
	beta.0.tmp=out.cond.2.full$beta.0.save[k]
	beta.tmp=out.cond.2.full$beta.save[,k]
  lam.save[,k]=exp(beta.0.tmp+X.sm%*%beta.tmp)
};cat("\n")
lam.mn=apply(lam.save,1,mean)
lam.u=apply(lam.save,1,quantile,.975)
lam.l=apply(lam.save,1,quantile,.025)

plot(s.sm,exp(beta.0+X.sm%*%beta),xlab="location",ylab=bquote(lambda),col=3,type="l")
polygon(c(s.sm,rev(s.sm)),c(lam.u,rev(lam.l)),col=rgb(0,0,0,.2),border=NA)
lines(s.sm,lam.mn,col=1,lwd=2)
