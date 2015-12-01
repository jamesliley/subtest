# Compute, plot and analyse single-SNP test statistics for between-subtype comparisons

##' Compute test statistic X1, posterior probability of Za and Zd in full model.
##'
##' @title X1
##' @param Z n x 2 matrix of Z scores; Z[,1]=Z_d, Z[,2]=Z_a
##' @param pars1 parameters of full model; in order, pi0, pi1, tau, sigma_1, sigma_2, rho
##' @return vector of values of X1
##' @export
##' @author James Liley
##' @examples
##' nn=100000
##' Z=abs(rbind(rmnorm(0.8*nn,varcov=diag(2)), rmnorm(0.15*nn,varcov=rbind(c(1,0),c(0,2^2))), rmnorm(0.05*nn,varcov=rbind(c(3^2,2),c(2,4^2)))));
##' pars=c(0.8,0.15,3,2,4,2)
##' X1=X1(Z,pars)
##' plotZ(Z,rlim=2); points(Z[which(X1>0.7),],col="red",pch=16)
X1=function(Z,pars1) getpx(Z,pars1)[,3] # error handlers covered in function getpx

##' Compute test statistic X2, difference in pseudolikelihood of each SNP between full and null model.
##'
##' @title X2
##' @param Z n x 2 matrix of Z scores; Z[,1]=Z_d, Z[,2]=Z_a
##' @param pars1 parameters of full model; in order, pi0, pi1, tau, sigma_1, sigma_2, rho
##' @param pars0 parameters of null model; in order, pi0, pi1, tau, sigma_1, sigma_2, rho
##' @param cond adjust for difference in likelihood of Z_a; 'conditional' pseudolikelihood difference
##' @return vector of values of X2
##' @export
##' @author James Liley
##' @examples
##' nn=100000
##' Z=abs(rbind(rmnorm(0.8*nn,varcov=diag(2)), rmnorm(0.15*nn,varcov=rbind(c(1,0),c(0,2^2))), rmnorm(0.05*nn,varcov=rbind(c(3^2,2),c(2,4^2)))));
##' pars1=c(0.8,0.15,3,2,4,2); pars0=c(0.8,0.1,3,4,1,0)
##' X=X2(Z,pars1,pars0,cond=TRUE); X[which(X<0)]=0; X[which(X>1)]=1
##' plotZ(Z,rlim=2,col=gray(X)); 
X2=function(Z,pars1,pars0,cond=TRUE) {

  # various error handlers
  if ((length(dim(Z))!=2)) stop("Z must be an n x 2 matrix, and 'weights' must be of length n")
  if (dim(Z)[2]!=2) stop("Z must be an n x 2 matrix")
  if (length(pars1)!=6) stop("Parameter 'pars1' must be a six-element vector containing values of (in order) pi0, pi1, tau, sigma1, sigma2, and rho")
  if (pars1[1]>=1 | pars1[2]>=1 | pars1[1]<=0 | pars1[2]<=0 | pars1[1]+pars1[2]>=1) stop("Values of pi0, pi1, and pi2 (pars1[1],pars1[2],1-pars1[1]-pars1[2]) must all be between 0 and 1")
  if (min(pars1[3:5])<=0 | pars1[6]<0) stop("The value of rho (pars1[6]) must be nonnegative and values of tau, sigma1, and sigma2 (pars1[3:5]) must be positive")
  if (pars1[6]>pars1[3]*pars1[5]) stop("The covariance matrix (sigma1^2 rho \\ rho sigma2^2) must be positive definite (rho < sigma1*sigma2)")
  if (length(pars0)!=6) stop("Parameter 'pars0' must be a six-element vector containing values of (in order) pi0, pi1, tau, sigma1, sigma2, and rho")
  if (pars0[1]>=1 | pars0[2]>=1 | pars0[1]<=0 | pars0[2]<=0 | pars0[1]+pars0[2]>=1) stop("Values of pi0, pi1, and pi2 (pars0[1],pars0[2],1-pars0[1]-pars0[2]) must all be between 0 and 1")
  if (min(pars0[3:5])<=0 | pars0[6]<0) stop("The value of rho (pars0[6]) must be nonnegative and values of tau, sigma1, and sigma2 (pars0[3:5]) must be positive")
  if (pars0[6]>pars0[3]*pars0[5]) stop("The covariance matrix (sigma1^2 rho \\ rho sigma2^2) must be positive definite (rho < sigma1*sigma2)")

dif=log(plhood(Z,pars1,sumlog=FALSE))-log(plhood(Z,pars0,sumlog=FALSE))

if (cond) dif=dif-log(plhood1(Z[,2],pars1,sumlog=FALSE))+log(plhood1(Z[,2],pars0,sumlog=FALSE))
dif
}


##' Compute test statistic X3, weighted geometric mean of Z_a and Z_d
##'
##' Defined as X3=Z_d^(alpha) Z_a^(1-alpha). Alpha is chosen to prioritise SNPs with high Z_d values (ie separating subtypes) over those with high Z_a (associated with phenotype. 
##'
##' Parameter alpha may be specified directly, or a 'good' value may be chosen from parameters. The value is chosen, given tau and sigma values, such that the points (tau,1) and (1,sigma) have the same value of X3. If no parameters are set, alpha defaults to 0.5. The parameter alpha is defined as a global variable at the end of the function.
##' 
##' @title X3
##' @param Z n x 2 matrix of Z scores; Z[,1]=Z_d, Z[,2]=Z_a
##' @param pars set of parameters of full model
##' @param tau value of tau used to compute alpha. Overrides parameter pars if set
##' @param sigma value of sigma used to compute alpha. Overrides parameter pars if set
##' @param alpha value of alpha set directly. Overrides parameters tau, sigma and pars if set. 
##' @return vector of values of X3 and set alpha as a global variable
##' @export
##' @author James Liley
##' @examples
##' nn=100000
##' Z=abs(rbind(rmnorm(0.8*nn,varcov=diag(2)), rmnorm(0.15*nn,varcov=rbind(c(1,0),c(0,2^2))), rmnorm(0.05*nn,varcov=rbind(c(3^2,2),c(2,4^2)))));
##' pars1=c(0.8,0.15,3,2,4,2); pars0=c(0.8,0.1,3,4,1,0)
##' X=X3(Z,pars1,pars0,cond=TRUE); X[which(X<0)]=0; X[which(X>1)]=1
##' plotZ(Z,rlim=2,col=gray(X)); 
X3=function(Z,pars=NULL,tau=NULL,sigma=NULL,alpha=NULL) {

  # various error handlers
if (dim(Z)[2]!=2) stop("Z must be an n x 2 matrix")

if (!is.null(pars)) {
  if (length(pars)!=6) stop("Parameter 'pars' must be a six-element vector containing values of (in order) pi0, pi1, tau, sigma1, sigma2, and rho")
  if (pars[1]>=1 | pars[2]>=1 | pars[1]<=0 | pars[2]<=0 | pars[1]+pars[2]>=1) stop("Values of pi0, pi1, and pi2 (pars[1],pars[2],1-pars[1]-pars[2]) must all be between 0 and 1")
  if (min(pars[3:5])<=0 | pars[6]<0) stop("The value of rho (pars[6]) must be nonnegative and values of tau, sigma1, and sigma2 (pars[3:5]) must be positive")
  if (pars[6]>pars[3]*pars[5]) stop("The covariance matrix (sigma1^2 rho \\ rho sigma2^2) must be positive definite (rho < sigma1*sigma2)")
}
if (!is.null(alpha) && (alpha>1 | alpha < 0)) stop ("Parameter alpha must be in (0,1)")
if ((is.null(tau) & !is.null(sigma)) | (!is.null(tau) & is.null(sigma))) stop("Either both or neither of tau and sigma should be set")

if (is.null(alpha)) {
 if (is.null(tau) & is.null(sigma)) {
  if (!is.null(pars)) {
   smax=max(pars[4:5]); tau=pars[3]
   alpha=log(smax)/(log(tau)+log(smax))
  } else alpha = 0.5 # if no parameters are set
 } else alpha=log(sigma)/(log(tau)+log(sigma))
}

alpha <<- alpha # define alpha as global variable
(Z[,1]^alpha)*(Z[,2]^(1-alpha))
}



##' Compute test statistic X4, conditional false discovery rate of Z_d given Z_a
##'
##' This test statistic is against a different null hypothesis than X1-X3; namely that the SNP of interest does not differentiate subgroups (ie, no requirement for high |Z_a|). 
##' 
##' This test statistic does account for Z_a, however, by implicitly adapting for correlation between Z_a and Z_d, effectively reducing the threshold for Z_d association for SNPs with high Z_a score, if there is evidence that both tend to be high concurrently.
##'
##' Note that the procedure of declaring non-null all SNPs with cFDR < alpha does not limit the false discovery rate of such SNPs to alpha(ie, the procedure is not analagous to the Benjamini-Hochberg FDR controlling procedure). A bound on the overall FDR can be obtained using the c2a function (below.)
##' @title X4
##' @param Z n x 2 matrix of Z scores; Z[,1]=Z_d, Z[,2]=Z_a
##' @param sub option to only calculate cFDR at a subset of Z scores; cFDR is computationally intensive to calculate.
##' @return vector of values of X4
##' @export
##' @author James Liley
##' @examples
##' nn=100000
##' Z=abs(rbind(rmnorm(0.8*nn,varcov=diag(2)), rmnorm(0.15*nn,varcov=rbind(c(1,0),c(0,2^2))), rmnorm(0.05*nn,varcov=rbind(c(3^2,2),c(2,4^2)))));
##' X=X4(Z,sub=which(Z[,1]^2 + Z[,2]^2 > 6))
##' plotZ(Z,rlim=2); points(Z[which(X<0.001),],col="red",pch=16)
X4=function(Z,sub=which(Z[,1]^2 + Z[,2]^2 > 6)) cfdr(Z,sub=sub)


##' Given some set of parameters in pars, compute the posterior likelihood of membership of each group for a set of SNP Z scores
##'
##' @title getpx
##' @param Z n x 2 matrix of Z scores; Z[,1]=Z_d, Z[,2]=Z_a
##' @param pars parameters of full model; in order, pi0, pi1, tau, sigma_1, sigma_2, rho
##' @param i get posterior probability of this group; if null, get posterior probability of all groups.
##' @return n x 3 matrix; element [i,j] is probability of SNP i being in group j.
##' @export
##' @author James Liley
##' @examples
##' nn=100000
##' Z=abs(rbind(rmnorm(0.8*nn,varcov=diag(2)), rmnorm(0.15*nn,varcov=rbind(c(1,0),c(0,2^2))), rmnorm(0.05*nn,varcov=rbind(c(3^2,2),c(2,4^2)))));
##' pars=c(0.8,0.15,3,2,4,2)
##' gp=getpx(Z,pars)
##' plotZ(Z,rlim=2,col=rgb(gp))
getpx=function(Z,pars) {

  # various error handlers
  if (dim(Z)[2]!=2) stop("Z must be an n x 2 matrix")
  if (length(pars)!=6) stop("Parameter 'pars' must be a six-element vector containing values of (in order) pi0, pi1, tau, sigma1, sigma2, and rho")
  if (pars[1]>=1 | pars[2]>=1 | pars[1]<=0 | pars[2]<=0 | pars[1]+pars[2]>=1) stop("Values of pi0, pi1, and pi2 (pars[1],pars[2],1-pars[1]-pars[2]) must all be between 0 and 1")
  if (min(pars[3:5])<=0 | pars[6]<0) stop("The value of rho (pars[6]) must be nonnegative and values of tau, sigma1, and sigma2 (pars[3:5]) must be positive")
  if (pars[6]>pars[3]*pars[5]) stop("The covariance matrix (sigma1^2 rho \\ rho sigma2^2) must be positive definite (rho < sigma1*sigma2)")
  
  pars=as.numeric(pars)
  ## likelihood for each group; faster having separate functions
  lhood1 =function(scale) 4*scale*exp(-rowSums(Z^2)/2)/(2*3.14159265)
  
  lhood2 =function(scale,sigma) 4*scale*exp(-(Z[,1]^2 + (Z[,2]/sigma)^2)/2)/(2*3.14159265*sigma)
  
  lhood3 <- function(scale,sigma,rho) {
    vc=rbind(c(sigma[1]^2,rho),c(rho,sigma[2]^2)); sv=solve(vc)
    vc2=rbind(c(sigma[1]^2,-rho),c(-rho,sigma[2]^2)); sv2=solve(vc2)
    2*scale*( exp(-rowSums((Z %*% sv) * Z)/2)/(2*3.14159265*sqrt(det(vc)))
              + exp(-rowSums((Z %*% sv2) * Z)/2)/(2*3.14159265*sqrt(det(vc2))))
  }

  px=cbind(lhood1(pars[1]),lhood2(pars[2],pars[4]),lhood3(1-pars[1]-pars[2],c(pars[3],pars[5]),pars[6]))
  px=px/rowSums(px)
  px
}




##' ###### Part of this description will go in the description for X4 and part of it in the description 
##' 
##' Compute the conditional false discovery rate (cFDR) at a set of paired Z scores.
##' 
##' Given two positively-valued random variables Z_i and Z_j, for which the distribution of Z_j is known and we have a hypothesis H on the distribution of Z_i, we define functions F_i and F_j as F_i(z_i)=Pr(Z_i>z_i|H), F_j(z_j)=Pr(Z_j>z_j) and set P_i, P_j as the random variables F_i(Z_i), F_j(Z_j) (these are analagous to p-values). For two observations z_i and z_j with p_i=F_i(z_i), p_j=F_j(z_j) we then define the conditional false discovery rate, or cFDR, as
##'  
##' cFDR(p_i|p_j)=Pr(H|P_i<p_i,P_j<p_j)
##' 
##' If P_i and P_j are independent under H, then the cFDR can be conservatively estimated as
##' 
##' cFDR = Pr(H) Pr(P_i<p_i|P_j<p_j,H)/Pr(P_i<p_i|P_j<p_j) < p_i*ord{p_j*|p_j*<p_j}/ord{(p_i*,p_j*)|(p_i*<p_i,p_j*<p_j)}
##'
##' that is, p_i divided by the proportion of observed (p_i*,p_j*) with p_j* less than p_j which also satisfy p_i*<p_i
##' 
##' This is a useful way to capture pleiotropy between phenotypes. In our case, if H is the hypothesis that Z_d~N(0,1), then Z_a and Z_d are independent.
##' 
##' Note that the procedure of declaring non-null all SNPs with cFDR < alpha does not limit the false discovery rate of such SNPs to alpha(ie, the procedure is not analagous to the Benjamini-Hochberg FDR controlling procedure). A bound on the overall FDR can be obtained using the c2a function (below.)
##' 
##' @title cfdr
##' @param Z an n x 2 vector of Z scores. Z_i is Z[,1], Z_j is Z[,2]
##' @param sub optional parameter; only compute cFDR at this subset of values.
##' @param #### Other stuff for shared control adjustment
##' @return vector of values of cFDR
##' @export
##' @author James Liley
##' @examples
##' nn=100000
##' Z=abs(rbind(rmnorm(0.8*nn,varcov=diag(2)), rmnorm(0.15*nn,varcov=rbind(c(1,0),c(0,2^2))), rmnorm(0.05*nn,varcov=rbind(c(3^2,2),c(2,4^2)))));
##' cfdr=cfdr(Z,sub=which(Z[,1]^2 + Z[,2]^2 > 6))
##' plotZ(Z,rlim=2); points(Z[which(cfdr<0.001),],col="red",pch=16)

cfdr=function(Z,sub=1:dim(Z)[1]) {

if (dim(Z)[2]!=2) stop("Z must be an n x 2 matrix")
if (!all(sub %in% 1:dim(Z)[1])) stop(paste0("Parameter sub must be a subset of 1:",dim(Z)[1]))
if (min(Z)<0) warning("Z contains negative values; two-sided p values used")
  
  p=2*pnorm(-Z[,1]); pc=2*pnorm(-Z[,2])
  cf=rep(1,length(p))
  for (i in 1:length(sub)) {
    ww=which(pc<pc[sub[i]]); 
    qq=(1+length(which(p[ww]<p[sub[i]])))/(1+length(ww))
    cf[sub[i]]=p[sub[i]]/qq
  }
  cf[which(cf>1)]=1
  cf
}




##' Compute the probability of X3 exceeding some threshold under a null hypothesis H0
##' 
##' H0 is defined as Z=(Z_d,Z_a) having a mixture distribution of N(0,I2) with probability pi0, N(0,(1,0 // 0,sigma)) with probability pi1 and N(0, (tau,0 // 0,1)) with probability pi2. Parameters pi1, pi2, tau and sigma can be derived from the fitted null model or set.
##' 
##' The p value is dependent on the value of alpha used to compute X3. If this is not specified in the function call, it is recalculated according to the same rules as the function X3
##' 
##' The probability is computed using a numerical integral over the (+/+) quadrant and the range and resolution of the integral can be set.
##' 
##' @title px3
##' @param X3c vector of X3 cutoffs at which to calculate probability
##' @param alpha value of alpha used in computing X3; recalculated if NULL
##' @param pars parameters in null model
##' @param pi two-element vector c(pi1,pi2); overrides pars if set. If pars and pi are null, defaults to (1/3, 1/3)
##' @param tau value in null distribution of Z; overrides pars if set. If pars and tau are null, defaults to 1.
##' @param sigma value in null distribution of Z; overrides pars if set. If pars and sigma are null, defaults to 1.
##' @param xmax compute integral over [0,xmax] x [0,xmax] as approximation to [0,inf] x [0,inf]
##' @param res compute integral at gridpoints with this spacing.
##' @return vector of p-values corresponding to X3c
##' @author James Liley
##' @export
##' @return list of probabilities
##' @examples
##' nn=100000
##' Z=abs(rbind(rmnorm(0.8*nn,varcov=diag(2)), rmnorm(0.1*nn,varcov=rbind(c(1,0),c(0,2^2))), rmnorm(0.1*nn,varcov=rbind(c(3^2,0),c(0,1))))); pars=c(0.8,0.1,2,3,1,0)
##' X=X3(Z,pars=pars)
##' subX=order(runif(nn))[1:500]
##' pX=px3(X[subX],alpha=alpha,pars=pars)
##' plot((1:length(pX))/(1+length(pX)),sort(pX),xlab="Quantile in U[0,1]",ylab="Observed quantile",main="Q-Q plot")

px3=function(X3c,alpha=NULL,pars=NULL,pi=NULL,tau=NULL,sigma=NULL,xmax=12,res=0.05) {
  
if (!is.null(pars)) {
  if (length(pars)!=6) stop("Parameter 'pars' must be a six-element vector containing values of (in order) pi0, pi1, tau, sigma1, sigma2, and rho")
  if (pars[1]>=1 | pars[2]>=1 | pars[1]<=0 | pars[2]<=0 | pars[1]+pars[2]>=1) stop("Values of pi0, pi1, and pi2 (pars[1],pars[2],1-pars[1]-pars[2]) must all be between 0 and 1")
  if (min(pars[3:5])<=0 | pars[6]<0) stop("The value of rho (pars[6]) must be nonnegative and values of tau, sigma1, and sigma2 (pars[3:5]) must be positive")
  if (pars[6]>pars[3]*pars[5]) stop("The covariance matrix (sigma1^2 rho \\ rho sigma2^2) must be positive definite (rho < sigma1*sigma2)")
}
if (!is.null(alpha) && (alpha>1 | alpha < 0)) stop ("Parameter alpha must be in (0,1)")
if (!(all(is.null(tau) & is.null(sigma) & is.null(pi)) | (all(!is.null(tau) & !is.null(sigma) & !is.null(tau))))) stop("Either all or none of pi, tau and sigma should be set")


# set alpha
if (is.null(alpha)) {
  if (is.null(tau) & is.null(sigma)) {
    if (!is.null(pars)) {
      smax=max(pars[4:5]); tau=pars[3]
      alpha=log(smax)/(log(tau)+log(smax))
    } else alpha = 0.5 # if no parameters are set
  } else alpha=log(sigma)/(log(tau)+log(sigma))
}

# set tau, sigma, pi1, pi2
if (is.null(pi) & is.null(sigma) & is.null(pi)) {
  if (!is.null(pars)) {
    pi=c(pars[2],1-pars[1]-pars[2])
    tau=pars[3]
    sigma=pars[4]
  } else {
    pi=c(1/3,1/3)
    tau=1
    sigma=1
  }
}


xm=seq(0,xmax,res); N=length(xm)
xmat=matrix(xm,N,N); ymat=t(xmat); Z=cbind(as.vector(xmat),as.vector(ymat))
lx=plhood(Z,c(1-pi[1]-pi[2],pi[1],tau,sigma,1,0),C=0,sumlog=FALSE)
# lx1=matrix(lhd(Z,y0d,C=0,sumlog=FALSE),N,N); image(log(lx1))
x3x=(Z[,1]^alpha) * (Z[,2]^(1-alpha)); lx=lx[order(-x3x)]
x3c=round((N^2) *(1-(ecdf(x3x)(X3c))))
cumsum(lx)[x3c]/sum(lx)

}





##' Compute an upper bound on the false discovery rate amongst SNPs with X4<X4c
##' 
##' Bound is based on estimating the area of the unit square satisfying X4<alpha. It is typically conservative.
##' 
##' Computation requires parametrisation of the distribution of Z_a, assumed to have a distribution of N(0,1) with probability pi0 and N(0,sigma^2) with probability 1-pi0. Values of pi0 and sigma can be obtained from the fitted parameters of the null model, or specified directly.
##' 
##' The probability is computed using a numerical integral over the (+/+) quadrant and the range and resolution of the integral can be set.
##' @title c2a
##' @param Z n x 2 matrix of Z scores; Z[,1]=Z_d, Z[,2]=Z_a
##' @param X4c vector of values of cFDR at which to compute overall FDR
##' @param pars0 parameters of null model; needed for determining P(Z_a<z_a)
##' @param pi0 proportion of SNPs not associated with phenotype; overrides pars0 if set
##' @param sigma standard deviation of observed Z_a scores in SNPs associated with the phenotype; overrides pars0 if set.
##' @param xmax compute integral over [0,xmax] x [0,xmax] as approximation to [0,inf] x [0,inf]
##' @param res compute integral at gridpoints with this spacing.
##' @param vector of FDR values corresponding to X4c
##' @author James Liley
##' @export
##' @return list of FDR values
##' @examples
##' nn=100000
##' Z=abs(rbind(rmnorm(0.8*nn,varcov=diag(2)), rmnorm(0.15*nn,varcov=rbind(c(1,0),c(0,2^2))), rmnorm(0.05*nn,varcov=rbind(c(3^2,2),c(2,4^2)))));
##' X=X4(Z,sub=which(Z[,1]^2 + Z[,2]^2 > 6))
##' Xm=which(X<0.05); Xsub=Xm[order(runif(length(Xm)))[1:100]] # sample of cFDR values 
##' true_fdr=rep(0,100); for (i in 1:100) true_fdr[i]=length(which(X[1:(0.95*nn)] <= X[Xsub[i]]))/length(which(X<=X[Xsub[i]])) # true FDR values (empirical)
##' fdr=c2a(Z,X[Xsub],pars0=c(0.8,0.15,3,2,4,2)) # estimated FDR using area method
##' plot(true_fdr,fdr,xlab="True FDR",ylab="Estimated",col="red"); points(true_fdr,X[Xsub],col="blue"); abline(0,1); legend(0.1*max(true_fdr),0.7*max(fdr),c("Area method", "cFDR"),col=c("red", "blue"),pch=c(1,1)) # cFDR underestimates true FDR; area method gives good approximation.
c2a=function(Z,X4c,pars0=NULL,pi0=NULL,sigma=NULL,xmax=12,ymax=12,res=0.01) {

if (!is.null(pars0)) {
  if (length(pars0)!=6) stop("Parameter 'pars0' must be a six-element vector containing values of (in order) pi0, pi1, tau, sigma1, sigma2, and rho")
  if (pars0[1]>=1 | pars0[2]>=1 | pars0[1]<=0 | pars0[2]<=0 | pars0[1]+pars0[2]>=1) stop("Values of pi0, pi1, and pi2 (pars0[1],pars0[2],1-pars0[1]-pars0[2]) must all be between 0 and 1")
  if (min(pars0[3:5])<=0 | pars0[6]<0) stop("The value of rho (pars0[6]) must be nonnegative and values of tau, sigma1, and sigma2 (pars0[3:5]) must be positive")
  if (pars0[6]>pars0[3]*pars0[5]) stop("The covariance matrix (sigma1^2 rho \\ rho sigma2^2) must be positive definite (rho < sigma1*sigma2)")
}
if (!(all(is.null(pi0) & is.null(sigma)) | (all(!is.null(pi0) & !is.null(sigma))))) stop("Either all or none of pi0 and sigma should be set")
  
if (is.null(pi0) & is.null(sigma)) {
  if (!is.null(pars0)) {
    pi0=pars0[1]+pars0[2]
    sigma=pars0[4]
  } else {
    pi0=0.8 # default values
    sigma=2
  }
}

gx=seq(0,xmax,res); gy=seq(0,ymax,res); gridx=matrix(gx,length(gx),length(gy)); gridy=t(matrix(gy,length(gy),length(gx))); 
Zx=cbind(as.numeric(gridx),as.numeric(gridy))

ww=which(Z[,1]^2 + Z[,2]^2 > 3.5^2)
qn=matrix(0,length(gx),length(gy)) # qn[i,j] is set to the number of elements of Z with Z[,1]>gx[i] and Z[,2]>gx[j]
e1=round(ecdf(gx)(Z[ww,1])*length(gx)); e2=round(ecdf(gy)(Z[ww,2])*length(gy)); nx=length(gx); ny=length(gy)
for (ii in 1:length(ww)) {
  qn[1:e1[ii],1:e2[ii]]=qn[1:e1[ii],1:e2[ii]]+1 # block out lower-left rectangle for every point
}
pj=2*pnorm(-Zx[,1]) # p value
qd=round(dim(Z)[1]*(1-ecdf(Z[,2])(Zx[,2]))) # denominator of observed quantile
cf=pj*(1+qd)/(1+as.vector(qn)); cf[which(Zx[,1]^2 + Zx[,2]^2 <3.6^2)]=1; cf=matrix(cf,length(gx),length(gy))

kk=4*((pi0*dmnorm(Zx,varcov=diag(2)))+ ((1-pi0)*dmnorm(Zx,varcov=rbind(c(1,0),c(0,sigma^2))))) # probability distribution of null SNPs

ycut=round(ecdf(cf)(X4c)*length(cf)) # number of cf values less than each value of X4c
null_L=cumsum(kk[order(cf)])[ycut]/sum(kk) # integrate pdf of null SNPs over region with cfdr less than X4c. Expected proportion of null SNPs lying in region with cfdr < cut.

# xij is defined such that xij[i,j]=probability of null SNP falling in region x>gx[i],y>gy[j]
xij=t(apply(apply(matrix(kk,length(gx),length(gy))[length(gx):1,length(gy):1], 2, cumsum), 1, cumsum))[length(gx):1,length(gy):1]
maxm=  cummax(as.numeric(xij)[order(cf)])[ycut]/sum(kk) # maxm[j] is the maximum expected number of null SNPs in a region (x>xx,y>yy) with cfdr(xx,yy)<X4c[j]
  
X4c*null_L/maxm # final FDR is bounded above by the ratio of expected number of SNPs in L divided by the expected number of SNPs in M

}
