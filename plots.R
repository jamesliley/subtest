##' Computes how well the parameters gamma and kappa are fitted at various numbers of simulations, and how these affect the computation of p-values. Returns an object of class 'sim_N_test' which can be plotted.
##'
##' @title Nsim_stat
##' @param lr matrix of simulation results
##' @param NN vector of integer; test performance of simulation fit for n simulations, where n takes the values in NN
##' @param NB number of bootstrap samples to take for each sample size in NN
##' @param p_at compute p values at each value in p_at
##' @param verbose print values of NN[i] as they complete.
##' @param quant_at for each size in NN and each parameter, record these quantiles of the bootstrapped parameter
##' @return list of two objects, stat and quant; stat[,ii,i] is, in order, kappa, gamma, and p-values given gamma and kappa of each value in p_at for the i'th bootstrap sample for bootstrap sample size NN[ii]; quant[k,ii,j] is the quant_at[j]'th quantile of boostrapped test statistic k amongst bootstrap samples of size NN[ii]
##' @author James Liley
##' @export
##' @examples
##' data("sim_mat.RData") 
##' stat=Nsim_stat(sim_mat)
##' summary(stat)
##' plot(stat)

Nsim_stat=function(lr,NN=c(50,100,500,1000,1500),NB=500,p_at=c(2,3,5,10),quant_at=pnorm(-2:2),verbose=TRUE) {

# Error handlers
if (dim(lr)[2]!=18) stop("Parameter lr must be a matrix of simulation results; each row should have the form [PL under H1] [PL under H0] [Correcting factor] [Number of iterations for H1] [Number of iterations of H0] [seed] [parameters under H1] [parameters under H0]")
if (!is.vector(NN)) stop("Parameter NN must be a vector of values")
if (dim(lr)[1]<max(NN)) warning("Largest value of NN is larger than total number of simulations; bootstrap samples of this size will be unreliable")


stat=array(0,dim=c(2+length(p_at),length(NN),NB))

for (ii in 1:length(NN)) {

  n=NN[ii]

  for (i in 1:NB) {
    boot=sample(1:dim(lr)[1],n,replace=TRUE)
    lr0=lr[boot,]
    S0=sim_analysis(lr0)
    
    stat[1,ii,i]=S0$kappa # kappa
    stat[2,ii,i]=S0$gamma # gamma
 
    for (j in 1:length(p_at)) stat[2+j,ii,i]=p_mixchi(p_at[j],S=S0)
  } 
  if (verbose) print(NN[ii])
}

quant=array(0,dim=c(dim(stat)[1:2],length(quant_at)))
for (i in 1:dim(stat)[1]) {
  for (j in 1:dim(stat)[2]) {
    quant[i,j,]=quantile(stat[i,j,],quant_at)
  }
}

yy=list(stat=stat,quant=quant,p_at=p_at,quant_at=quant_at,B_size=NN)

class(yy)="sim_N_test"
yy
}


##' Print method for class sim_N_test
##' @title print.Nsim_stat
##' @param yy object of class sim_N_test, output from Nsim_stat
##' @return for each bootstrap sample size, prints quantiles of each parameter and then all data.
##' @author James Liley
##' @export
print.sim_N_test=function(yy) {
cat("Number of bootstrap samples per run: ",dim(yy$stat)[3],"\n\n")
for (ii in 1:length(yy$B_size)) {
  cat("Bootstrap sample size: ",yy$B_size[ii],"\n")
  ds=data.frame(t(yy$quant[,ii,]))
  rownames(ds)=round(yy$quant_at,digits=3); colnames(ds)=c("kappa","gamma",paste0("Pr(PLR>",yy$p_at,")"))
  print(ds)
  cat("\n")
}
for (ii in 1:length(yy$B_size)) {
  cat("Bootstrap sample size: ",yy$B_size[ii],"\n")
  ds=data.frame(t(yy$stat[,ii,]))
  colnames(ds)=c("kappa","gamma",paste0("Pr(PLR>",yy$p_at,")"))
  print(ds)
  cat("\n")
}
}

##' Summary method for class sim_N_test
##' @title summary.Nsim_stat
##' @param yy object of class sim_N_test, output from Nsim_stat
##' @return for each bootstrap sample size, prints quantiles of each parameter
##' @author James Liley
##' @export
summary.sim_N_test=function(yy) {
  cat("Number of bootstrap samples per run: ",dim(yy$stat)[3],"\n\n")
  for (ii in 1:length(yy$B_size)) {
    cat("Bootstrap sample size: ",yy$B_size[ii],"\n")
    ds=data.frame(t(yy$quant[,ii,]))
    rownames(ds)=round(yy$quant_at,digits=3); colnames(ds)=c("kappa","gamma",paste0("Pr(PLR>",yy$p_at,")"))
    print(ds)
    cat("\n")
  }
}


##' Plot method for class sim_N_test
##' 
##' Draws plot of values and quantiles of each test statistic for each bootstrap sample size
##' 
##' @title summary.Nsim_stat
##' @param yy object of class sim_N_test, output from Nsim_stat
##' @param all combine all plots in one
##' @author James Liley
##' @export
plot.sim_N_test=function(yy,all=FALSE) {
  
if (all) par(mfrow=c(1,dim(yy$stat)[1])) else par(mfrow=c(1,1))
  
xx=as.numeric(t(matrix(1:length(yy$B_size),length(yy$B_size),dim(yy$stat)[3])))
  
pc_lf=c("p","p",rep("f",length(yy$p_at))) # plot percentage or fold change
names=c(expression(kappa),expression(gamma),paste0("Pr(PLR > ",yy$p_at,")"))

for (i in 1:(2+length(yy$p_at))) {
  
  if (pc_lf[i]=="p") {
    q2=100*((yy$quant[i,,]/median(yy$quant[i,,]))-1); 
    s2=100*((yy$stat[i,,]/median(yy$stat[i,,]))-1)
  } else{
    q2=log10(yy$quant[i,,]/median(yy$quant[i,,])); 
    s2=log10(yy$stat[i,,]/median(yy$stat[i,,]));
  }
  if (all) ylab0="" else {
    if (pc_lf[i]=="p") ylab0="Percentage error" else ylab0="Fold error"
  }
  plot(0,0,type="n",xlim=c(1,length(yy$B_size)),ylim=range(q2),main=names[i],xaxt="n",xlab=if (all) "" else "Size of bootstrap samples",ylab=ylab0,las=1,yaxp=c(round(range(q2)),4))
  for (ii in 1:dim(yy$quant)[3]) {
    lines(1:dim(yy$quant)[2],q2[,ii],col=c("blue","red","white","red","blue")[ii])
  }
  abline(h=0,lty=2)
  legend(1.2,0,signif(median(yy$stat[i,,]),digits=2),yjust=0.5,bg="white",box.col="white")
  
  y=as.numeric(t(s2)); 
  points(xx,y,pch=".",col="darkgray")
  
  if (!all) {
    axis(1,at=1:length(yy$B_size),labels=yy$B_size)
    if (i<dim(yy$stat)[1]) { cat ("Press [enter] to continue");    line <- readline() }
  }
}
if (all) {
  axis(1,at=1:length(yy$B_size),labels=yy$B_size)
  mtext("Error", side=2, at=17,line=3,cex=0.8)
  par(mfrow=c(1,1))
}
}  

