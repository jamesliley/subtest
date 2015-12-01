# Functions for analysing random subtypes

##' For a genotype matrix, simulate a random subtyping, compute corresponding Za and Zd scores, and fit model.
##'
##' @title rand_gen
##' @param X a snpMatrix object, QC'd, or a list of objects of type snpMatrix covering different SNPs (ie, one per chromosome). X is assumed to only contain cases as rows, as it somewhat reduces efficiency to split the matrix every time.
##' @param Z_a Z_a scores (case vs control); these will be the same for each random subgroup (as the case set is the same) so it is inefficient to recalculate each time.
##' @param pars_init1 an nx6 array of starting parameter values for fitting full model. The E-M algorithm is run from each of them in turn.
##' @param pars_init0 an nx6 array of starting parameter values for fitting null model. The E-M algorithm is run from each of them in turn.
##' @param weights SNP weights to adjust for LD; output from LDAK procedure
##' @param Ca matrix of covariates for computing Z_a; Ca[i,] gives the covariates for sample i corresponding to X[i,]
##' @param Cd matrix of covariates for computing Z_d; Cd[i,] gives the covariates for sample i corresponding to X[i,]. CD[i,] is ignored unless Ya[i]>0. 
##' @param n0 optionally, only simulate subtypes such that the smaller subtype has size n0. If n0 is null, use random subgroup sizes with the smaller subtype has between 10% and 50% of total cases.
##' @param Yd optionally choose random subtyping by randomly sorting 'true' Yd. This is the best option for quantitative traits Yd. Supercedes n0 if set. 
##' @param file set to a directory where output of the simulation will be saved. The output file is small but the simulation is likely to be repeated multiple times. If file is null, print output of the simulation to console.
##' @param seed random seed for computing random subtypes. Can be used to reconstruct simulation. By default, chosen by time.
##' @param ... other arguments passed to fit.3g
##' @export
##' @author James Liley
##' @examples
##' # pending

rand_gen=function(X,Z_a,pars_init1=c(0.8,0.15,2,1,2,1),pars_init0=c(0.8,0.15,2,2,1,0),weights=rep(1,length(Z_a)),Ca=NULL,Cd=NULL,n0=NULL,Yd=NULL, file=NULL,seed=NULL,...) {

# Error handlers
#
# Errors in the form of X and Cd will be caught  by function zd_scores

if (dim(X)[2]!=length(Z_a)|dim(X)[2]!=length(weights)) stop("Parameter Z_a must correspond to the columns of X. Length of Z_a must be the same as the number of columns of X and the length of parameter 'weights'")
if (length(pars_init1)!=6| !is(pars_init1,"numeric")| max(pars_init1[1:2])>1| pars_init1[1]+pars_init1[2]>1| min(pars_init1[1:5])<= 0 | pars_init1[6]==0 | pars_init1[6]>pars_init1[3]*pars_init1[5])
 stop("Parameter pars_init1 must be a six-element vector containing elements (pi0,pi1,tau,sigma_1,sigma_2,rho). The first five elements must be strictly positive and the sixth nonnegative. Parameters pi0 and pi1 must be less than 1 and sum to less than 1. Parameter rho must be less than tau*sigma_2.")
if (length(pars_init0)!=6| !is(pars_init0,"numeric")| max(pars_init0[1:2])>1| pars_init0[1]+pars_init0[2]>1| min(pars_init0[1:5])<= 0 | pars_init0[6]==0 | pars_init0[6]>pars_init0[3]*pars_init0[5])
 stop("Parameter pars_init0 must be a six-element vector containing elements (pi0,pi1,tau,sigma_1,sigma_2,rho). The first five elements must be strictly positive and the sixth nonnegative. Parameters pi0 and pi1 must be less than 1 and sum to less than 1. Parameter rho must be less than tau*sigma_2.")
if (n0>dim(X)[1]) stop("Parameter n0 must be less than the number of rows of X")
if (!(is.null(Yd)||length(Yd)!=dim(X)[1])) stop("If specified, parameter Yd must have the same length as the number of rows of X")

if (is.null(seed)) {
  options(digits.secs=6);
  seed=as.numeric(substr(Sys.time(),21,26))
}  

set.seed(seed) # set random seed

if (!is.null(Yd)) {
  Yp=Yd[order(runif(dim(X)[1]))] # set response variable by random parameter 
} else {
  if (is.null(n0)) n0=round(runif(1,0.1,0.5)*dim(X)[1])
  Yp=c(rep(0,n0),rep(1,dim(X)[1]-n0))[order(runif(dim(X)[1]))]
}

# if (length(unique(Yp))>2) fam="Gaussian" else fam="Binomial"

zd=zd_scores(X,Yp,Cd=Cd,signed=FALSE)
#if (is(X,"SnpMatrix") | is(X,"XSnpMatrix")) zd=z_scores1(X,Yp,Cd=Cd,signed=FALSE,fam=fam) else {
#  zd=c()
#  for (i in 1:length(X)) zd=c(zd,z_scores1(X[[i]],Yp,Cd=Cd,signed=FALSE,fam=fam))
#}
z_ad=abs(cbind(zd,Z_a))

# clean data
ww=which(!is.na(z_ad[,1]+z_ad[,2]))
z_ad=z_ad[ww,];
weights=weights[ww]


parsf1=matrix(0,dim(pars_init1)[1],6) # fitted parameters
lhd1=rep(0,dim(pars_init1)[1]) # likelihood of Z_d, Z_a|pars
lha1=rep(0,dim(pars_init1)[1]) # likelihood of Z_a|pars
n1=rep(0,dim(pars_init1)[1]) # number of iterations of EM algorithm for convergence
for (i in 1:dim(pars_init1)[1]) {
  yy1=fit.3g(z_ad,pars=pars_init1[i,],weights=weights,...)
  parsf1[i,]=yy1$pars
  lhd1[i]=yy1$logl
  lha1[i]=yy1$logl_a
  n1[i]=dim(yy1$hist)[1]
}

parsf0=matrix(0,dim(pars_init0)[1],6)
lhd0=rep(0,dim(pars_init0)[1])
lha0=rep(0,dim(pars_init0)[1])
n0=rep(0,dim(pars_init0)[1])

for (i in 1:dim(pars_init0)[1]) {
  yy0=fit.3g(z_ad,pars=pars_init0[i,],fit_null=TRUE,weights=weights,...)
  parsf0[i,]=yy0$pars
  lhd0[i]=yy0$logl
  lha0[i]=yy0$logl_a
  n0[i]=dim(yy0$hist)[1]
}


wx0=which.max(lhd0)
pars0=parsf0[wx0,]
lh0=lhd0[wx0]
la0=lha0[wx0]
n0x=n0[wx0]

wx1=which.max(lhd1)
pars1=parsf1[wx1,]
lh1=lhd1[wx1]
la1=lha1[wx1]
n1x=n1[wx1]

vec=c(lh0,lh1,la1-la0,n0x,n1x,seed,pars0,pars1)

if (!is.null(file)) write(vec,paste0(file,"/sim",as.character(seed),".txt"),append=FALSE,ncolumns=18) else print(vec)

}


##' For a given directory or matrix of outputs from function sim_gen, computes a corresponding set of start points for the E-M algorithm on subsequent simulations.
##'
##' @title rand_start
##' @param sims either a matrix of outputs from sim_gen, a directory, or a list of files.
##' @param n generate n start points for simulation.
##' @param mindist generate start points for simulation at least this far apart (distance measured by d_par). Overrides n if set.
##' @param nsim use this many simulations to determine cluster points. Reading a large number of files (as is the case when parameter 'sims' is a list) can be slow.
##' @param maxit if the E-M algorithm uses the maximum allowed number of iterations, it may not have converged. Set this parameter to restrict to only simulations taking fewer than this number of iterations.
##' @export
##' @author James Liley
##' @examples
##' data(sim_mat)
##' sim_start(sim_mat)

rand_start=function(sims,n=2,mindist=NULL,nsim=100,maxit=1e4) {

# error handlers
if (!(is.matrix(sims)|is.data.frame(sims)|is.character(sims))) stop("Parameter sims must be a matrix, a list of files, or a directory containing simulation output results")
if (is.character(sims)) {
  if (length(sims)==1 && !file.exists(sims)) stop(paste0("Directory ",sims," not found"))
  if (length(sims)>1) for (i in 1:length(sims)) if (!file.exists(sims[i])) stop(paste0("File ",sims[i]," not found"))
}  

if (is.matrix(sims) | is.data.frame(sims)) {
  nsim=min(nsim,dim(sims)[1])
  pars=sims[(1:dim(sims)[1])[order(runif(dim(sims)[[1]]))[1:nsim]],] # random selection of nsim rows
}
if (is.character(sims) & length(sims)>1) {
  nsim=min(nsim,length(sims))
  sims_shuffle=sims[order(runif(length(sims)))] # random ordering
  pars=c()  
  for (i in 1:nsim) pars=rbind(pars,read.table(sims_shuffle[i]))
}
if (is.character(sims) & (length(sims)==1)) {
  files=list.files(sims,pattern="sim*")
  nsim=min(nsim,length(files))
  sims_shuffle=files[order(runif(length(files)))]
  pars=c()
  for (i in 1:nsim) pars=rbind(pars,read.table(paste0(sims,"/",sims_shuffle[i])))
}

sims=pars

if (dim(sims)[2]!=18) stop("Each file or row of simulation matrix must have a single line with 18 elements; in order, the PL under H1, the PL under H0, the correction for Za, the number of iterations of the EM algorithm in the full and null models, the random seed, the paramaters of the full model, and the parameters of the null model.")

wnc=which(sims[,4]<maxit & sims[,5]<maxit)
if (length(wnc)<dim(sims)[1]) message(paste0(dim(sims)[1]-length(wnc)," simulations removed due to non-convergence"))
sims=sims[wnc,]

if (dim(sims)[1]<n) stop(paste0("Too few valid simulations for ",n," clusters"))
   
pars1=pars[,13:18]; pars0=pars[,7:12]

if (!is.null(mindist)) cc1=cutree(hclust(d_par(pars1)),h=mindist) else cc1=cutree(hclust(d_par(pars1)),k=n)
if (!is.null(mindist)) cc0=cutree(hclust(d_par(pars0)),h=mindist) else cc0=cutree(hclust(d_par(pars0)),k=n)

xpars0=c()
for (i in 1:length(unique(cc0))) xpars0=rbind(xpars0,colMeans(pars0[which(cc0==i),]))

xpars1=c()
for (i in 1:length(unique(cc1))) xpars1=rbind(xpars1,colMeans(pars1[which(cc1==i),]))

return(list(pars1=xpars1,pars0=xpars0))

}






##' Given outputs of fitted parameters and likelihoods from random subtypes, analyse outcomes in order to compute p-values by extrapolation or interpolation.
##'
##' @title rand_analysis
##' @param sims either a matrix of outputs from sim_gen, a directory, or a list of files.
##' @param maxit if the E-M algorithm uses the maximum allowed number of iterations, it may not have converged. Set this parameter to restrict to only simulations taking fewer than this number of iterations.
##' @export
##' @author James Liley
##' @examples
##' data(sim_mat)
##' S=sim_analysis(sim_mat)
##' plot(S)
rand_analysis=function(sims,maxit=1e4) {

# error handlers
if (!(is.matrix(sims)|is.data.frame(sims)|is.character(sims))) stop("Parameter sims must be a matrix, a list of files, or a directory containing simulation output results")
if (is.character(sims)) {
  if (length(sims)==1 && !file.exists(sims)) stop(paste0("Directory ",sims," not found"))
  if (length(sims)>1) for (i in 1:length(sims)) if (!file.exists(sims[i])) stop(paste0("File ",sims[i]," not found"))
}  
  
if (is.matrix(sims) |is.data.frame(sims)) X=sims
if (is.character(sims) & length(sims)>1) {
  X=c()  
  for (i in 1:length(sims)) X=rbind(X,read.table(sims[i]))
}
if (is.character(sims) & (length(sims)==1)) {
  files=list.files(sims,pattern="sim*")
  X=c()
  for (i in 1:length(files)) X=rbind(X,read.table(paste0(sims,"/",files[i])))
}

if (dim(X)[2]!=18) stop("Each file or row of simulation matrix must have a single line with 18 elements; in order, the PL under H1, the PL under H0, the correction for Za, the number of iterations of the EM algorithm in the full and null models, the random seed, the paramaters of the full model, and the parameters of the null model.")

wnc=which(X[,4]<maxit & X[,5]<maxit)
if (length(wnc)<dim(X)[1]) message(paste0(dim(X)[1]-length(wnc)," simulations removed due to non-convergence"))
X=X[wnc,]


N=dim(X)[1]

plr=X[,2]-X[,1]-X[,3]; plr[which(plr<0)]=0
tau=X[,15]; s1=X[,16]; s2=X[,17]; rho=X[,18]
  

yx=mixchipars(plr) # get gamma and kappa

yy=list(plr=plr,tau=tau,s1=s1,s2=s2,rho=rho,kappa=yx[2],gamma=yx[1])

class(yy)="sim_output"
yy
}

##' Fit maximum-likelihood values of gamma and kappa to a set of PLR values
##' 
##' Assumes that PLR has a mixture distribution with some weight p0 at 0, a weight 'kappa' for a distribution gamma*chi_1^2, and a weight 1-kappa-p0 for a distribution gamma*chi_2^2.
##' 
##' @title mixchipars
##' @param plr set of plr values
##' @param init values of gamma, kappa at which to start numerical fitting algorithm
##' @return two element vector; gamma, kappa.
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


##' Print method for class sim_output
##' @param yy object of class sim_output; generally output from sim_analysis
##' @export
##' @author James Liley
print.sim_output=function(yy) {
  cat("Fitted values\n")
  cat("kappa",yy$kappa,"\n")
  cat("gamma",yy$gamma,"\n")
  cat("\n")
  print(data.frame(PLR=yy$plr,tau=yy$tau,sigma_1=yy$s1,sigma_2=yy$s2,rho=yy$rho))
}


##' Summary method for class sim_output
##' @param yy object of class sim_output; generally output from sim_analysis
##' @author James Liley
summary.sim_output = function(yy) {
  cat("Fitted values\n")
  cat("kappa",yy$kappa,"\n")
  cat("gamma",yy$gamma,"\n")
  cat("\n")
  
  plr=yy$plr; tau=yy$tau; s1=yy$s1; s2=yy$s2; rho=yy$rho
  min_s=rep(0,length(s1)); min_s[which(s1<s2)]=s1[which(s1<s2)]; min_s[which(s1>s2)]=s2[which(s1>s2)]; 
  max_s=rep(0,length(s1)); max_s[which(s1<s2)]=s2[which(s1<s2)]; max_s[which(s1>s2)]=s1[which(s1>s2)]; 
  
  cat("PLR\n")
  cat("Min",min(plr),"\n")
  cat("Quantile 0.25",quantile(plr,0.25),"\n")
  cat("Median",median(plr),"\n")
  cat("Quantile 0.75",quantile(plr,0.75),"\n")
  cat("Max",max(plr),"\n")
  cat("\n")
  
  cat("Tau\n")
  cat("Min: ",min(tau),"\n")
  cat("Max: ",max(tau),"\n")
  cat("Mean: ",mean(tau),"\n")
  cat("St. dev: ",sd(tau),"\n")
  cat("Proportion <1: ",length(which(tau<1))/length(tau),"\n")
  cat("\n")

  cat("Sigma_1/Sigma_2\n")
  cat("Prop. sig_1 > s_2: ",length(which(s1>s2))/length(s1),"\n")
  cat("Mean min(s1,s2): ",mean(min_s),"\n")
  cat("St. dev. min(s1,s2): ",mean(min_s),"\n")
  cat("Mean max(s1,s2): ",mean(max_s),"\n")
  cat("St. dev. max(s1,s2): ",mean(max_s),"\n")
  cat("Max s1,s2: ",max(max_s))
  cat("\n")
  
  cat("Rho")
  cat("Prop. zero: ",length(which(rho<0.01))/length(rho),"\n")
  cat("Median: ",median(rho),"\n")
  cat("Max: ",max(rho),"\n")
  cat("\n")
  
}

##' Plot method for class sim_output. Plots the values of PLR in ascending order along with mixture chi-squared distribution, and a Q-Q plot comparing the two. Sets global variables X,Y where (X,Y) are the points on the Q-Q plot.
##' @param yy object of class sim_output; generally output from sim_analysis
##' @param QQ set to TRUE to only draw QQ plot, FALSE to only draw chi-squared plot, or NULL to draw both.
##' @param conf set to TRUE to draw a 95% confidence interval on the Q-Q plot. Confidence intervals are constructed by repeatedly randomly sampling from a mixture chi-squared distribution parametrised by yy$gamma and yy$kappa. For the random sample, the best-fit values of gamma and kappa (gamma' and kappa') are recovered, and the random sample is compared to the fitted mixture-chi square parametrised by gamma' and kappa'. Finally, standard errors are estimated empirically at each quantile. Computing limits is slow, so if this option is TRUE, global variables lb and ub are set, where (X, lb), (X,ub) are the co-ordinates of the lower and upper confidence limits. 
##' @param q95 set to TRUE to draw lines on the plot corresponding to empirical type 1 error if generating a significance cutoff from mixture-chi squared corresponding to p<0.05.
##' @param ... additional parameters passed to plot
##' @author James Liley
plot.sim_output=function(yy,QQ=NULL, conf=FALSE,q95=FALSE,...) {


if ((!is.null(QQ) & !is.logical(QQ)) | !is.logical(conf)) stop ("Parameters QQ and conf must be a logical value or NULL")

plr=yy$plr
N=length(plr)
NN=10000; n0=round(NN*length(which(plr==0))/length(plr)); n1=round(yy$kappa*NN); n2=NN-n0-n1
ya=yy$gamma*sort(c(rep(0,n0),qchisq((1:n1)/(1+n1),df=1),qchisq((1:n2)/(1+n2),df=2))); ya0=quantile(ya,(1:N)/(1+N))

if (is.null(QQ) || !QQ) {
 plot((1:N)/(1+N),sort(plr),xlab="Relative rank",ylab="PLR",...); lines((1:NN)/(1+NN),ya,lwd=2,col="red")
 legend(0.2,0.8*max(plr),c("Observed",expression(paste("Mixture ",chi^2))),pch=c(1,45),col=c("black","red"),bty="n")
}

if (is.null(QQ)) {
 cat ("Press [enter] to continue")
 line <- readline()
}

if (is.null(QQ)|| QQ) {
 plot(ya0,sort(plr),type="n",xlab=expression(paste("Expected quantile (mixture-",chi^2,")")),ylab="Observed PLR",xaxs="i",yaxs="i",...); abline(0,1,col="red",lwd=3)
 
 if (conf) { # draw confidence limits
   if (!exists("lb")|!exists("ub")||!is.numeric(lb)||!is.numeric(ub)||length(lb)!=length(plr)||length(ub)!=length(plr)) { # ie, if lb and ub are not already set.
   N=500 # perform this many simulations to generate limits
   r_q=matrix(0,length(ya0),N) # random sample from mixture-chi2 parametrised by yy$gamma and yy$kappa
    
   nn0=length(which(plr==0)); nn1=round(length(plr)*yy$kappa); nn2=length(plr)-nn0-nn1 
   for (i in 1:N) {
     ys=sort(c(rep(0,nn0),rchisq(nn1,df=1),rchisq(nn2,df=2)))*yy$gamma
     gk=mixchipars(ys);
     rg=gk[1];rk=gk[2] # maximum-likelihood estimates of gamma and kappa on random dataset
     nx1=round(length(plr)*rk); nx2=length(plr)-nn0-nx1 ;
     yx=sort(c(rep(0,nn0),rchisq(nx1,df=1),rchisq(nx2,df=2)))*rg
     r_q[,i]=approx(yx,ys,ya0)$y
   }
   
   lb=rep(0,length(ya0)); ub=rep(0,length(ya0))
   for (i in 1:length(ya0)) {
     lb[i]=quantile(r_q[i,],0.05,na.rm=TRUE)
     ub[i]=quantile(r_q[i,],0.95,na.rm=TRUE)
   }
   lb <<- lb; ub <<- ub # set as global variables
   }
   polygon(c(ya0,max(ya0),ya0[length(ya0):1]),c(lb,max(ya0),ub[length(ub):1]),col="grey",border=NA)
 }
 
 points(ya0,sort(plr),cex=0.5); abline(0,1,col="red",lwd=3)
 X<<- ya0; Y<<- sort(plr)
 if (q95) {
   qx=quantile(ya0,0.95); qy=sort(plr)[length(which(ya0<qx))]
   abline(v=qx,lty=2,col="blue")
   abline(h=qy,lty=2,col="blue")
   text(0.5*(qx+max(ya0)),0.8*qy,paste0("Obs. type 1 err: ",signif(length(which(plr>qx))/length(plr),3)))
 }
 legend(0.2,0.8*max(plr),c("Observed",expression(paste("Mixture-",chi^2))),pch=c(1,45),col=c("black","red"),bty="n") 

}
}

##' Computes the p-value of an observed PLR according to either empirical quantile in simulations, a mixture chi-square distribution, or the empirical bound, depending on size.
##' 
##' @title p_value
##' @param x list of plr values at which to compute p value
##' @param S object of class sim_output; output from sim_analysis
##' @author James Liley
##' @export
##' @examples
##'  gamma=0.7; kappa=0.3; N=1000; 
##'  n1=round(kappa*N); n2=N-n1; Xs=gamma*c(qchisq((1:n1)/(1+n1),df=1),qchisq((1:n2)/(1+n2),df=2))
##'  
##'  pp=p_mixchi(Xs,gamma=gamma,kappa=kappa)
##'  plot((1:length(pp))/(1+length(pp)),sort(pp),xlab="Quantile in U(0,1)",ylab="P-value"); abline(0,1,col="red")
##' p_mixchi(10,gamma,kappa)
p_value=function(x,S=NULL) {

if (is.null(S)) p_mixchi(x) else { # if S is null, p_mixchi will throw a warning message because gamma and kappa are not set
out=rep(1,length(x))
out[which(x<=0)]=1

if (length(S$plr)<200) cut1=quantile(S$plr,0.9) else cut1=sort(S$plr)[length(S$plr)-20] # below this cutoff, use empirical quantile in PLR
cut2=8*quantile(S$plr,0.9) # below this cutoff, use mixture chi-squared; above, use MLE-based upper bound

w1=which(x<cut1) # at this range, it is OK to use empirical quantile as p-value
if (length(w1)>0) out[w1]=1-ecdf(S$plr)(x[w1])


w2=which(x>= cut1 & x<cut2) # roughly, our simulations indicate that the chi-squared approximation is valid over this range.
if (length(w2)>0) out[w2]=p_mixchi(x[w2], S=S)

w3=which(x>= cut2) # in this range, use robust p-value bound
if (length(w3)>0) {
  out[w3]=pbind(x[w3],S=S)
  out[w3][which(out[w3]>p_mixchi(cut2,S=S))]=p_mixchi(cut2,S=S) # ensure p value is a non-increasing function of PLR.
}
}

return(out)
}


##' Computes a p-value assuming the PLR follows a mixture chi-squared distribution.
##' 
##' @title p_mixchi
##' @param x list of plr values at which to compute p value
##' @param S object of class sim_output; output from sim_analysis
##' @param gamma scaling factor for chi-squared. Overrides sim if set.
##' @param kappa mixing proportion of chi_1^2. Overrides sim if set
##' @param p0 proportion of plr values which are negative. Overrides sim if set. If gamma and kappa are set, default is 0.
##' @return list of p-values
##' @author James Liley
##' @export
##' @examples
##'  gamma=0.7; kappa=0.3; N=1000; 
##'  n1=round(kappa*N); n2=N-n1; Xs=gamma*c(qchisq((1:n1)/(1+n1),df=1),qchisq((1:n2)/(1+n2),df=2))
##'  
##'  pp=p_mixchi(Xs,gamma=gamma,kappa=kappa)
##'  plot((1:length(pp))/(1+length(pp)),sort(pp),xlab="Quantile in U(0,1)",ylab="P-value"); abline(0,1,col="red")
##' p_mixchi(10,gamma,kappa)
p_mixchi=function(x,S=NULL,gamma=NULL, kappa=NULL,p0=NULL) {

if ((is.null(gamma) & !is.null(kappa)) | (!is.null(gamma) & is.null(kappa))) stop("Either both or neither of parameters gamma and kappa must be set")
if (!is.null(S) && !is(S,"sim_output")) stop ("Parameter S must be of class sim_output")
if (!is.null(gamma) & !is.null(kappa) & is.null(p0)) p0=0

if (is.null(gamma) & is.null(kappa)) {
  if (!is.null(S)) {
    gamma=S$gamma
    kappa=S$kappa
    p0=length(which(S$plr==0))/length(S$plr); if (!is.finite(p0)) p0=0
  } else {
    gamma=0.5 # default value
    kappa=0.5
    warning("Parameters gamma and kappa not set. P value will be unreliable")
  }
}

(1-p0)*((kappa*pchisq(x/gamma,df=1,lower.tail=FALSE)) + ((1-kappa)*pchisq(x/gamma,df=2,lower.tail=FALSE)))

}



##' Computes a bound on the p-value for a subtype comparison, relatively robust to LD, by principally using fitted values rather than overall likelihood. 
##' 
##' The bound is very conservative for most observed pseudo-likelihood ratios. However, at very high likelihood ratios, the mixture-chi squared approximation appears to break down, meaning that p-values generated using this algorithm could result in underestimated false-discovery rates. We recommend that this function is used to estimate p-values if the observed PLR is more than twice the maximum obtained in simulated data, or for PLR greater than the maximum simulated value if the observed distribution of PLR statistics appears to poorly match the mixture-chi squared distribtion at high PLR values.
##'
##' @title pbind
##' @param x list of plr values at which to compute p value
##' @param S object of class sim_out, containing PLR, tau, s1, s2, rho. Overrides other parameters if set.
##' @param plr list of plr values from simulation
##' @param tau list of fitted values of tau corresponding to values of plr
##' @param s1 list of fitted values of sigma_1 corresponding to values of plr
##' @param s2 list of fitted values of sigma_2 corresponding to values of plr
##' @param rho list of fitted values of rho corresponding to values of plr
##' @param res resolution of numerical integral 
##' @author James Liley
##' @export
##' @examples lr=read.table("~/Subtypes/AI_data/simulation/LR");
##' plr=lr[,2]-lr[,1]-lr[,3]; plr[which(plr<0)]=0; tau=lr[,34]; s1=lr[,36]; s2=lr[,37]; rho=lr[,38]
##' pbind(c(0.5,1,1.5,2,3,5,10,20,30),plr,tau,s1,s2,rho)
pbind=function(x,S=NULL,plr=NULL,tau=NULL,s1=NULL,s2=NULL,rho=NULL,res=1000) {
if (is.null(S) & (is.null(plr)|is.null(tau)|is.null(s1)|is.null(s2)|is.null(rho))) stop("Either parameter S or all of parameters plr, tau, s1, s2, and rho must be set")
if (!(class(S) %in% c("NULL","sim_output"))) stop("Parameter S must be of class sim_output")

if (!is.null(S)) {
  plr=S$plr
  tau=S$tau
  s1=S$s1
  s2=S$s2
  rho=S$rho
}

x_init=x

if (max(x_init)<quantile(plr,0.9)) ecdf(plr)(x_init) else { # if p value is high, estimate directly
  
  w=which(x_init<quantile(plr,0.9)); 
  if (length(w)>0) {
    y=x_init[w]; 
    p_big=ecdf(plr)(y)
    x=x_init[-w]
  } else x=x_init
  
  
  xi=length(which(tau>1))/length(tau) # proportion of simulations with tau>1
  nu=length(which(tau>1 & s2>s1))/length(which(tau>1))
  
  # Case where tau<1
  q1=which(tau<1)
  
  adj0=mean(sort(plr[q1])[1:2]); yy=plr[q1]-adj0; 
  
  kappa=length(which(tau<1 & rho==0))/length(which(tau<1)) # proportion of simulations with tau<1 and rho=0
  N=10000; n0=round(N*kappa); n1=N-n0; ya=sort(c(qchisq((1:n0)/(1+n0),df=1),qchisq((1:n1)/(1+n1),df=2)))
  beta=lm(sort(yy)~quantile(ya,(1:length(yy))/(1+length(yy))))$coefficients # PLR ~ beta[2]*(kappa*chi[1]^2 + (1-kappa)*chi[2]^2); y  =yabeta
  ya=ya*beta[2] + beta[1]
    
  # if (length(which(ya>x-adj0))>0) p_t1=length(which(ya>(x-adj0)))/length(ya) else p_t1=pchisq((x-adj0-beta[1])/beta[2],df=2,lower.tail=FALSE)
  
  p_t1=1-ecdf(ya)(x-adj0)
  qc2=which(x>max(ya)+adj0); p_t1[qc2]=pchisq((x[qc2]-adj0-beta[1])/beta[2],df=2,lower.tail=FALSE)
  
  
  # Case where s1>s2, tau>1
  p1=which(s1>s2 & tau>1)
  adj1=mean(sort(plr[p1])[1:2]); yy=plr[p1]-adj1; xx=rho[p1]^4; rc=lm(yy~xx)$coefficients # plr = rc*rho^4
  
  x1=rho[p1][which(rho[p1]>0)]^2
  sr1=sqrt(mean(x1^2)) # standard deviation of normal describing positive part of rho distribution
  sn0=length(which(rho[p1]>0))/length(p1) # proportion of instances of rho which are > 0
  
  rho2_cut=sqrt((x-adj1-rc[1])/rc[2]) # rho^2 has to reach this cutoff for plr to be >x
  p_s1=2*sn0*(1-pnorm(rho2_cut,sd=sr1)) # approximate probability of rho^2 reaching this value
  
  
  
  # Case where s1<s2, tau>1
  p2=which(s2>s1 & tau>1)
  t2=tau[p2]
  r2=(rho[p2]/(tau[p2]*s2[p2]))^2
  
  st2=sqrt(mean((t2-1)^2))
  sr2=sqrt(mean(r2[which(r2>0)]^2))
  pr20=length(which(r2>0))/length(r2) # probability that r^2 > 0
  
  fft=(t2^2) - log(t2^2)
  ffr=log(1-r2)
  adj2=mean(sort(plr[p2])[1:2]); xx=(fft-ffr-1); yy=plr[p2]-adj2
  mm=quantile(yy/xx,0.99) #max((yy/xx)[which(yy>quantile(yy,0.5))])
  
  xxcut=(x-adj2)/mm # broadly, xx must meet this cutoff for plr>x
  
  
  # To work out the bound on Pr(PLR>x|s2>s1,tau>1), we evaluate a double integral over the region (tau,r^2) in [0,inf] x [0,1]
  
  # Auxiliary function: solve e^(x-c)=x by Newton-Raphson method
  tau_solve=function(c,y1,tol=1e-5) {
    y0=0; ii=0; while (abs(y1-y0) > tol & ii<20) {
      ii=ii+1; y0=y1; y1=y0-((exp(y0-c) - y0)/(exp(y0-c)-1)) 
    }; if (ii<20) y1 else 2*c
  }
  
  tmax=sqrt(tau_solve(max(xxcut)+1,(max(xxcut)+1)*1.3)) # If tau>tmax, then F(tau,r^2)>x whatever the value of r^2
  ptmax=2*(1-pnorm(tmax,mean=1,sd=st2)) # probability of tau exceeding this value
  
  
  # integrate over this many boxes in x and y directions
  xr=1+((tmax-1)*(1:res)/(1+res)); pxr=2*dnorm(xr,mean=1,sd=st2)
  yr=(1:res)/(1+res); pyr=pr20*2*dnorm(yr,sd=sr2)
  
  probs=outer(pxr,pyr) # value of pdf of tau,r^2 at grid points
  pass=outer((xr^2)-log(xr^2),log(1-yr),function(x,y) x-y-1) # value of F(tau,r^2)
  
  probs1=as.vector(probs)
  pass1=as.vector(pass)
  
  passcuts=round((1-ecdf(pass1)(xxcut))*length(pass1))
  
  probs2=probs1[order(-pass1)]
  
  p_s2=ptmax + (cumsum(probs2)[passcuts]*(tmax-1)/(res^2))
  
  p_small=as.numeric(
    (xi*( (nu*p_s1) + ((1-nu)*p_s2) )) + ((1-xi)*p_t1)
  )
  
  if (length(x_init)==1) return(p_small) else {
    p_all=rep(0,length(x_init))
    if (length(w)>0) {
      p_all[w]=p_big
      p_all[-w]=p_small
    } else p_all=p_small
    return(p_all)
  }
  
}

}
  