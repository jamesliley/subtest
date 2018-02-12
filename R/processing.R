# Functions to process a (QC'ed) genotype matrix for three-group test.


##' Obtain a set of starting parameters for E-M algorithm, given Z_a and Z_d scores and LDAK weights. Because of the potential for local maxima in the likeilhood landscape, it is important to start the algorithm at several points. Since the E-M algorithm is computationally intensive, it is useful to start as close as possible to the actual maxima, and at as few points as possible. This function seeks a small number of 'promising' start points with high likelihood, sufficiently far from each other.
##' 
##' @title pars_start
##' @param Z an n x 2 array; Z[i,1], Z[i,2] are the Z_d and Z_a scores respectively for the ith SNP
##' @param weights SNP weights to adjust for LD; output from LDAK procedure
##' @param H hypothesis, 0 or 1
##' @param n1 begin with n1 random parameter sets from pars_rand, chosen according to prior distributions of parameters (see documentation for pars_rand)
##' @param n2 trim initial list to this many well-separated sets of parameters; each of these parameter sets undergoes nx steps of the EM algorithm
##' @param n3 finally take n3 well-separated sets of parameters for entry into the final E-M algorithm
##' @param nx use this many iterations of the E-M algorithm on each of the n2 sets of paramaters
##' @param C scaling factor for adjustment
##' @param seed random seed for generating results. Use to regenerate.
##' @author James Liley
##' 
pars_start=function(Z,weights=rep(1,dim(Z)[1]),H=1,n1=100,n2=10,n3=5,C=1,nx=3,seed=NULL) {

  # various error handlers
  if ((length(dim(Z))!=2) | (dim(Z)[1]!=length(weights))) stop("Z must be an n x 2 matrix, and 'weights' must be of length n")
  if (dim(Z)[2]!=2) stop("Z must be an n x 2 matrix")
  if (!((n1 >= n2) & (n2 >= n3))) stop("Parameters n1, n2, n3 must be successively smaller.")
  
  
pars1=pars_rand(n1,seed=seed)
ll=apply(pars1,1,function(x) plhood(Z,x,weights=weights,C=C))

dp=d_par(pars1)
hp=hclust(dp) # hiererarchical clustering of parameters
cp=cutree(hp,n2) # select n2 clusters of parameters

pars2=matrix(0,n2,6)
for (i in 1:n2) {
  sub=which(cp==i);
  if (length(sub)==1) pars2[i,]=pars1[sub,] else pars2[i,]=pars1[sub[which.max(ll[sub])],]
} # best set of parameters in each cluster

pars2a=pars2
ll3=rep(0,n2)
for (i in 1:n2) {
  yy=fit.3g(Z,pars2[i,], weights=weights,C=1,maxit=nx,accel=FALSE,fit_null=(H==0))
  pars2a[i,]=yy$pars
  ll3[i]=yy$logl
} # run each cluster through 10 steps of the algorithm

cp3=cutree(hclust(d_par(pars2a)),k=n3) # cluster again

pars3=matrix(0,n3,6)
for (i in 1:n3) {
  sub=which(cp3==i);
  if (length(sub)==1) pars3[i,]=pars2a[sub,] else pars3[i,]=pars2a[sub[which.max(ll3[sub])],]
} # best set of parameters in each cluster

pars3
}



##' Randomly generate a set of 'feasible' parameter sets describing distribution of Z_a and Z_d, used for initiating E-M algorithm.
##' pi2 is drawn such that log10(pi2) has a uniform distribution with max. value 5, min. value 1. Values of pi2>0.1 are not expected, and 10^-5 is a low enough prior to initialise the EM algorithm if the true value is smaller than this.
##' pi1 is drawn from either U(0.3,0.7) or such that log10(pi1)~U(1,3), the first with probability p_reg
##' sigma_1 and tau are either 1 or drawn such that sigma_1^2~U(1,10)
##' sigma_2 is either 1 or drawn such that sigma_1^2~U(1,10)
##' rho is either 0 or U(0.05,0.95)*sigma_2*tau
##' @name pars_rand
##' @param N number to generate 
##' @param H hypothesis; 0 or 1. Under H=0, sigma2 is forced to 1 and rho to 0.
##' @param p_reg probability that pi1 is U(0.3,0.7), sigma_1=1, sigma_2=1, tau=1, rho=0
##' @param seed random seed for generating results; use to regenerate.
##' @author James Liley
##' @export
##' @examples
##' px=pars_rand()
##' par(mfrow=c(2,3))
##' plot(density(px[,2],from=0,to=1),main="pi1") # distribution of pi1
##' plot(density(log10(1-px[,1]-px[,2]),from=-5,to=-1),main="log10(pi2)") # empirical distribution of log10(pi2)
##' plot(density(px[,3]),main="tau") # empirical density of tau
##' plot(density(px[,4]),main="sigma1") # empirical density of sigma1
##' plot(density(px[,5]),main="sigma2") # empirical density of sigma2
##' plot(density(px[,6],from=0),main="rho") # empirical density of rho
##' par(mfrow=c(1,1))
pars_rand = function(N=1000,H=1,p_reg=0.3,seed=NULL) {

if (is.null(seed)) {
  options(digits.secs=6);
  seed=as.numeric(substr(Sys.time(),21,26)) # choose from clock time
}
set.seed(seed)

  if (p_reg<=0 |p_reg >= 1) stop("Parameter p_reg must be in (0,1)")
  if (!(H %in% c(0,1))) stop("H must be 0 or 1")

  n0=round(p_reg*N); n1=N-n0
  p2=10^(-runif(N,1,5)) # max 0.1
  p1=c(runif(n0,0.3,0.7),10^(-runif(n1,1,3)))*(1-p2) # can be either approx. 0.5
  p0=1-p1-p2
  
  t=c(rep(1,n0),sqrt(runif(n1,0,10)))[order(runif(N))]
  s1=c(rep(1,n0),sqrt(runif(n1,0,10)))[order(runif(N))]
  if (H==1) {
    s2=c(rep(1,n0),sqrt(runif(n1,0,16)))[order(runif(N))]
    rho=c(rep(0,n0),runif(n1,0.05,0.95))[order(runif(N))]*s2*t
  } else {
    s2=rep(1,N)
    rho=rep(0,N)
  }
  cbind(p0,p1,t,s1,s2,rho)
}





##' Distance metric between sets of parameters. Computed as a Euclidean distance in a transformed parameter space. pi0-pi2 are transformed to their inverse logits, as a difference 0.001-0.1 means more than a difference 0.4-0.5. Tau, sigma1, and sigma2 are transformed to sqrt(|x-1|), as differences are most important for values near 1. Rho is transformed to x^2 /2, since differences are only important for large rho.
##' @title d_par
##' @param pars1 an n x 6 matrix of parameters
##' @param ... additional parameters passed to dist
##' 
##' @return an n x m matrix A such that A[i,j] is the distance between pars1[i,] and pars2[j,]
##' @author James Liley
##' @export
##' @examples
##' px=pars_start(N=5)
##' dx=d_par(px)
##' image(as.matrix(dx))
d_par=function(pars1,...) {
  ## inverse logit
  il=function(x) -log((1/x)-1)
  
  pars1=as.matrix(pars1);
  p1=cbind(il(pars1[,1:2]),il(1-pars1[,1]-pars1[,2]),sqrt(abs(pars1[,3:5]-1)),0.5 * (pars1[,6]^2))
  dist(p1,...)
}
