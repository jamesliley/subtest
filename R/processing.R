# Functions to process a (QC'ed) genotype matrix for three-group test.


##' Generates Z_a and Z_d scores from a QC'd genotype matrix
##'
##' @title z_scores
##' @param X a snpMatrix object, QC'd, or a list of objects of type snpMatrix covering different SNPs (ie, one per chromosome)
##' @param Ya indices of cases/controls; Ya[i] corresponds to X[i,]. 0 for controls, 1 for cases
##' @param Yd indices of subgroups; Yd[i] corresponds to X[i,]. NA for controls, 1 for subtype 1, 2 for subtype 2, or the value of a quantitative trait across cases.
##' @param Ca matrix of covariates for computing Z_a; Ca[i,] gives the covariates for sample i corresponding to X[i,]
##' @param Cd matrix of covariates for computing Z_d; Cd[i,] gives the covariates for sample i corresponding to X[i,]. CD[i,] is ignored unless Ya[i]>0. 
##' @param signed set to TRUE to return signed Z scores
##' @param control set to TRUE to perform genomic control on Z scores. We recommend strict genomic control (ie, so median(p)=0.5).
##' @export
##' @author James Liley
##' @examples
##' # pending - need to find a SnpMatrix which can be made publically available

z_scores=function(X,Ya,Yd,Ca=NULL,Cd=NULL,signed=TRUE,control=TRUE) {

if (is(X,"SnpMatrix")) X=list(X=X) # make handling easier
nm=length(X) # number of SnpMatrices

# various error handlers
if (!("snpStats" %in% rownames(installed.packages()))) stop("Requires package snpStats")
for (i in 1:nm) {
  if (!(is(X[[i]],"SnpMatrix") | is(X[[i]],"XSnpMatrix"))) stop("X must be a SnpMatrix/XSnpMatrix or list of objects of these classes (package snpStats)")
  if (!((dim(X[[i]])[1]==length(Ya)) & (dim(X[[i]])[1]==length(Yd)) & (is.null(Ca) || (dim(Ca)[1]==dim(X[[i]])[1])) & (is.null(Cd) || (dim(Cd)[1]==dim(X[[i]])[1])))) stop("Lengths of Ya and Yd, number of rows of Ca and Cd, and number of rows of each SnpMatrix must be the same")
}
if (!(is.null(Ca) | is.matrix(Ca) |is.data.frame(Ca) |is.null(Cd) | is.matrix(Cd) |is.data.frame(Cd))) stop("Parameters Ca and Cd must be matrices or data frames")
  
# Required scripts and packages
require(snpStats,quiet=TRUE)
  
names(Ya)=rownames(X[[1]]); if (!is.null(Ca)) rownames(Ca)=rownames(X[[1]]);
names(Yd)=rownames(X[[1]]); if (!is.null(Cd)) rownames(Cd)=rownames(X[[1]]);
  
z_ad=cbind(zd_scores(X,Yd,Cd,signed,control=control),za_scores(X,Ya,Ca,signed,control=control))
z_ad
}





##' Generates Z_a scores from a QC'd genotype matrix
##'
##' @title za_scores
##' @param X a snpMatrix object, QC'd, or a list of objects of type snpMatrix covering different SNPs (ie, one per chromosome)
##' @param Ya indices of cases/controls; Ya[i] corresponds to X[i,]. 0 for controls, 1 for cases
##' @param Ca matrix of covariates for computing Z_a; Ca[i,] gives the covariates for sample i corresponding to X[i,]
##' @param signed set to TRUE to return signed Z scores
##' @param control set to TRUE to perform genomic control on Z scores
##' @export
##' @author James Liley
##' @examples
##' # See examples for function z_scores
za_scores=function(X,Ya,Ca=NULL,signed=TRUE,control=TRUE) {

if (is(X,"SnpMatrix")) X=list(X=X) # make handling easier
nm=length(X) # number of SnpMatrices

# various error handlers
if (!("snpStats" %in% rownames(installed.packages()))) stop("Requires package snpStats")
for (i in 1:nm) {
  if (!(is(X[[i]],"SnpMatrix") | is(X[[i]],"XSnpMatrix"))) stop("X must be a SnpMatrix/XSnpMatrix or list of objects of these classes (package snpStats)")
  if (!((dim(X[[i]])[1]==length(Ya)) & (is.null(Ca) || (dim(Ca)[1]==dim(X[[i]])[1])))) stop("Length of Ya, number of rows of Ca, and number of rows of each SnpMatrix must be the same")
  }
if (!(is.null(Ca) | is.matrix(Ca) |is.data.frame(Ca))) stop("Parameter Ca must be a matrix or data frame")

names(Ya)=rownames(X[[1]]); if (!is.null(Ca)) rownames(Ca)=rownames(X[[1]]); # row name set; required for SnpStats

za=c()
for (i in 1:nm) za=c(za,z_scores1(X[[i]],Ya,Ca,signed))
if (control) {
  qa=qchisq(2*pnorm(-abs(za)),df=1,lower.tail=FALSE); 
  lambda_a=median(qa,na.rm=TRUE)/qchisq(0.5,df=1,lower.tail=FALSE)
  qa=qa/lambda_a
  za=-sign(za)*qnorm(pchisq(qa,df=1,lower.tail=FALSE))
}
za
}




##' Generates Z_d scores from a QC'd genotype matrix
##'
##' @title zd_scores
##' @param X a snpMatrix object, QC'd, or a list of objects of type snpMatrix covering different SNPs (ie, one per chromosome)
##' @param Yd indices of subgroups; Yd[i] corresponds to X[i,]. 1 for subtype 1, 2 for subtype 2, or the value of a quantitative trait across cases. MUST be NA for controls and cases for which subtype is unknown or undefined.
##' @param Cd matrix of covariates for computing Z_d; Cd[i,] gives the covariates for sample i corresponding to X[i,]. Cd[i,] is ignored if Yd[i]=NA
##' @param signed set to TRUE to return signed Z scores
##' @param control set to TRUE to perform genomic control on Z scores##' @export
##' @author James Liley
##' @examples
##' # See examples for function z_scores
zd_scores=function(X,Ya,Ca=NULL,signed=TRUE,control=TRUE) {

if (is(X,"SnpMatrix")) X=list(X=X) # make handling easier  
nm=length(X) # number of SnpMatrices

# various error handlers
if (!("snpStats" %in% rownames(installed.packages()))) stop("Requires package snpStats")
for (i in 1:nm) {
  if (!(is(X[[i]],"SnpMatrix") | is(X[[i]],"XSnpMatrix"))) stop("X must be a SnpMatrix/XSnpMatrix or list of objects of these classes (package snpStats)")
  if (!((dim(X[[i]])[1]==length(Ya)) & (is.null(Ca) || (dim(Ca)[1]==dim(X[[i]])[1])))) stop("Length of Yd, number of rows of Cd, and number of rows of each SnpMatrix must be the same")
}
if (!(is.null(Cd) | is.matrix(Cd) |is.data.frame(Cd))) stop("Parameter Cd must be a matrix or data frame")

names(Yd)=rownames(X[[1]]); if (!is.null(Cd)) rownames(Cd)=rownames(X[[1]]); # set row names; required for SnpStats

w=which(!is.na(Yd))
if (length(w)==length(Yd)) {
 sYd=Yd
 if (!is.null(Cd)) sCd=Cd else sCd=NULL 
} else { # Remove NA elements (controls)
 sYd=Yd[w]; 
 if (!is.null(Cd)) sCd=Cd[w,] else sCd=NULL # sub-indices and matrices only including cases
}

if (length(unique(sYd))>2) fam="Gaussian" else fam="Binomial" # use a gaussian family if Yd is interval rather than categorical


zd=c()
for (i in 1:length(X)) {
  if (length(w)==length(Yd)) sX=X[[i]] else sX=X[[i]][w,];
  zd=c(zd,z_scores1(sX,sYd,sCd,signed,fam))
}
if (control) {
  qa=qchisq(2*pnorm(-abs(za)),df=1,lower.tail=FALSE); 
  lambda_a=median(qa,na.rm=TRUE)/qchisq(0.5,df=1,lower.tail=FALSE)
  qa=qa/lambda_a
  za=-sign(za)*qnorm(pchisq(qa,df=1,lower.tail=FALSE))
}

zd
}

# Auxiliary function for single SnpMatrix
z_scores1=function(X,Yd,Cd=NULL,signed=TRUE,fam="Binomial") {
  if (!is.null(Cd)) {
    form=as.formula(paste("Yd~",paste(colnames(sCd),collapse=" + ")))
    xd=snp.rhs.tests(Yd~as.matrix(Cd[w,]),data=Cd[w,],snp.data=X,score=signed,fam=fam) 
    if (signed) sign=effect.sign(xd,simplify=TRUE) else sign=1
    zd=-qnorm(p.value(xd)*sign/2)
  } else {
    xd=single.snp.tests(Yd,snp.data=X,score=signed)
    if (signed) sign=effect.sign(xd) else sign=1
    zd=-qnorm(p.value(xd,df=1)/2)*sign
  }
  zd
}



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
pars_start=function(Z,weights=rep(1,dim(Z)[1]),H=1,n1=1000,n2=10,n3=3,C=1,nx=10,seed=NULL) {

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
