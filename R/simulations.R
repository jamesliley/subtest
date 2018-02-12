# Functions for analysing random subtypes


##' Computes the p-value of an observed PLR according to either empirical quantile in simulations, a mixture chi-square distribution, or the empirical bound, depending on size.
##' 
##' @title p_value
##' @param x list of plr values at which to compute p value
##' @param cplr list of cplr values from simulations from which to estimate null distribution
##' @author James Liley
##' @export
##' @examples
p_value=function(x,cplr) {

out=rep(1,length(x))
out[which(x<=0)]=1

p=x

for (i in 1:length(x))
  if (x[i]<quantile(cplr,0.9)) p[i]=ecdf(cplr)(x[i]) else {
    mm=mixchipars(cplr)
    gamma=mm[1]
    kappa=mm[2]
    p0=length(which(cplr==0))/length(cplr); if (!is.finite(p0)) p0=0
    p[i]=(1-p0)*((kappa*pchisq(x[i]/gamma,df=1,lower.tail=FALSE)) +
                   ((1-kappa)*pchisq(x[i]/gamma,df=2,lower.tail=FALSE)))


  }

p
}


##' Fit maximum-likelihood values of gamma and kappa to a set of PLR values
##' 
##' Assumes that PLR has a mixture distribution with some weight p0 at 0, a weight 'kappa' for a distribution gamma*chi_1^2, and a weight 1-kappa-p0 for a distribution gamma*chi_2^2.
##' 
##' @title mixchipars
##' @param plr set of plr values
##' @param init values of gamma, kappa at which to start numerical fitting algorithm
##' @return two element vector; gamma, kappa.
##' @export
##' @examples
##' gamma=0.7;kappa=0.4;p0=0.05; N=10000
##' px=sort(gamma*c(rep(0,p0*N),rchisq(kappa*N,df=1),rchisq((1-kappa)*N,df=2)))
##' mixchipars(px)
mixchipars=function(plr,init=c(1,0.5))  {
  plr1=plr[which(plr>0)]
  plr_lhd=function(gk) if (gk[1]<=0 |gk[2]<0|gk[2]>1) 1e50 else min(-sum(log((gk[2]*dchisq(plr1/gk[1],df=1)/gk[1]) + ((1-gk[2])*dchisq(plr1/gk[1],df=2)/gk[1]))),1e50) # Negative likelihood of observing a particular set of values of PLR given kappa (gk[2]) and gamma (gk[1]). Must not be infinite.
  
  yx=optim(init,plr_lhd)$par #yx=optim(c(1,0.5),plr_lhd,lower=c(1e-3,0),upper=c(100,1),method="L-BFGS-B")$par # optimise over a 
  yx[2]=yx[2]*length(plr1)/length(plr) # scale since plr1 is only positive values of plr.
  yx
}

