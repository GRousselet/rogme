#' Difference asymmetry function for two independent groups
#'
#' The difference asymmetry function provides perspective on the degree a distribution
#' is symmetric about zero, by quantifying the sum of q and 1-q quantiles of the
#' distribution of all pairwise differences between two independent groups x and y.
#'  If the distributions of x and y are identical, then the
#'  distribution of all pairwise differences is symmetric about zero.
#' If the distribution is symmetric the function should be approximately a horizontal line.
#' If in addition the median of the difference scores is zero, the horizontal line will
#' intersect the y-axis at zero.
#' Confidence intervals and p values are returned for each quantile sum.
#' The FWE is controlled via Hochberg's method, which is used to determine critical
#' p values based on the argument alpha.
#'
#' @section Note:
#' This function combines Rand Wilcox's qwmwhd & cbmhd R functions,
#' from Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/},
#' http://dornsife.usc.edu/labs/rwilcox/software/.
#' @references
#' Wilcox, R.R., Erceg-Hurn, D.M., Clark, F. & Carlson, M. (2014)
#' Comparing two independent groups via the lower and upper quantiles.
#' J Stat Comput Sim, 84, 1543-1551.
#' @seealso \code{\link{hd}}, \code{\link{diff_asym}}
asymhd <- function(){

}

qwmwhd<-function(x,y,q=seq(5,40,5)/100,xlab="Quantile",ylab="Sum of q and 1-q Quantiles",plotit=TRUE,alpha=.05,nboot=1000){
  #
  #  Plot that provides perspective on the degree a distribution is symmetric about zero.
  #  This function plots the sum of q and 1-q quantiles of the distribution of D=X-Y, X and Y independent.
  #  A 1-alpha confidence interval for the sum is indicated by a +
  #  If the distribution is symmetric
  #  the plot should be approximately a horizontal line.
  #
  #  FWE is controlled via Hochberg's method, which was used to determine critical
  #  p-values based on the argument
  #  alpha.
  #
  #  Can alter the quantiles compared via the argument
  #  q
  #  q must be less than .5
  #
  x=elimna(x)
  y=elimna(y)
  n1=length(x)
  n2=length(y)
  output=matrix(NA,ncol=8,nrow=length(q))
  dimnames(output)=list(NULL,c("quantile","Est.1","Est.2","SUM","ci.low","ci.up","p_crit","p-value"))
  for(i in 1:length(q)){
    test=cbmhd(x,y,q=q[i],plotit=FALSE,nboot=nboot)
    output[i,1]=q[i]
    output[i,2]=test$Est1
    output[i,3]=test$Est2
    output[i,4]=test$sum
    output[i,8]=test$p.value
    output[i,5]=test$ci[1]
    output[i,6]=test$ci[2]
  }
  temp=order(output[,8],decreasing=TRUE)
  zvec=alpha/c(1:length(q))
  output[temp,7]=zvec
  output <- data.frame(output)
  output$signif=rep("YES",nrow(output))
  for(i in 1:nrow(output)){
    if(output[temp[i],8]>output[temp[i],7])output$signif[temp[i]]="NO"
    if(output[temp[i],8]<=output[temp[i],7])break
  }
  if(plotit){
    plot(rep(q,3),c(output[,4],output[,5],output[,6]),type="n",xlab=xlab,ylab=ylab)
    points(q,output[,6],pch="+")
    points(q,output[,5],pch="+")
    points(q,output[,4],pch="*")
  }
  list(n=c(n1,n2),output=output)
}

cbmhd<-function(x,y,alpha=.05,q=.25,plotit=FALSE,pop=0,fr=.8,rval=15,xlab="",ylab="",nboot=600,SEED=TRUE){
  #
  #  Compute a confidence interval for the sum of the qth and (1-q)th quantiles
  #  of the distribution of D=X-Y, where X and Y are two independent random variables.
  #  The Harrell-Davis estimator is used
  #  If the distribution of X and Y are identical, then in particular the
  #  distribution of D=X-Y is symmetric about zero.
  #
  #  plotit=TRUE causes a plot of the difference scores to be created
  #  pop=0 adaptive kernel density estimate
  #  pop=1 results in the expected frequency curve.
  #  pop=2 kernel density estimate (Rosenblatt's shifted histogram)
  #  pop=3 boxplot
  #  pop=4 stem-and-leaf
  #  pop=5 histogram
  #
  if(SEED)set.seed(2)
  if(q>=.5)stop("q should be less than .5")
  if(q<=0)stop("q should be greater than 0")
  x<-x[!is.na(x)]
  y<-y[!is.na(y)]
  n1=length(x)
  n2=length(y)
  m<-outer(x,y,FUN="-")
  q2=1-q
  est1=hd(m,q)
  est2=hd(m,q2)
  data1<-matrix(sample(n1,size=n1*nboot,replace=TRUE),nrow=nboot)
  data2<-matrix(sample(n2,size=n2*nboot,replace=TRUE),nrow=nboot)
  bvec=NA
  for(i in 1:nboot){
    mb=outer(x[data1[i,]],y[data2[i,]],"-")
    bvec[i]=hd(mb,q)+hd(mb,q2)
  }
  p=mean(bvec>0)+.5*mean(bvec==0)
  p=2*min(c(p,1-p))
  sbv=sort(bvec)
  ilow<-round((alpha/2) * nboot)
  ihi<-nboot - ilow
  ilow<-ilow+1
  ci=sbv[ilow]
  ci[2]=sbv[ihi]
  if(plotit){
    if(pop==1 || pop==0){
      if(length(x)*length(y)>2500){
        print("Product of sample sizes exceeds 2500.")
        print("Execution time might be high when using pop=0 or 1")
        print("If this is case, might consider changing the argument pop")
        print("pop=2 might be better")
      }}
    MM=as.vector(m)
    if(pop==0)akerd(MM,xlab=xlab,ylab=ylab)
    if(pop==1)rdplot(MM,fr=fr,xlab=xlab,ylab=ylab)
    if(pop==2)kdplot(MM,rval=rval,xlab=xlab,ylab=ylab)
    if(pop==3)boxplot(MM)
    if(pop==4)stem(MM)
    if(pop==5)hist(MM,xlab=xlab)
    if(pop==6)skerd(MM)
  }
  list(q=q,Est1=est1,Est2=est2,sum=est1+est2,ci=ci,p.value=p)
}

# ==========================================================================

#' Difference asymmetry function for two dependent groups or a single distribution
#'
#' \code{diff_asym} computes a difference asymmetry function for two dependent
#' groups or a single distribution. The difference asymmetry function provides
#' perspective on the degree a distribution is symmetric about zero, by
#' quantifying the sum of q and 1-q quantiles. If the distribution is symmetric
#' the function should be approximately a horizontal line. If in addition the
#' median of the difference scores is zero, the horizontal line will intersect
#' the y-axis at zero. Confidence intervals and p values are returned for each
#' quantile sum. The FWE is controlled via Hochberg's method, which is used to
#' determine critical p values based on the argument alpha.
#'
#' @param x A a vector of pairwise differences.
#' @param q Quantiles to estimate - default = 0.05:0.05:.4 - must be <.5.
#' @param nboot Number of bootstrap samples - default = 1000.
#' @param alpha Expected long-run type I error rate - default = 0.05.
#'
#' @return
#' saved in structure diff_asym_res with fields:
#' - q = estimated quantiles
#'- hd_q = Harrell-Davis estimate of quantile q
#'- hd_1minusq = Harrell-Davis estimate of quantile 1-q
#'- qsum = sum of hd_q and hd_1minusq - the main quantity of interest
#'- qsumci = matrix quantiles x low/high bounds of the confidence
#'          intervals of the quantile sums
#'- pval_crit = critical p value, corrected for multiple comparisons
#'- pval = p value
#'
#' @section Note:
#' This function combines Rand Wilcox's difQpci and Dqdif R functions,
#' from Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/},
#' http://dornsife.usc.edu/labs/rwilcox/software/.
#' @references
#' Wilcox, R.R. & Erceg-Hurn, D.M. (2012)
#' Comparing two dependent groups via quantiles.
#' J Appl Stat, 39, 2655-2664.
#' #' @seealso \code{\link{hd}}, \code{\link{diffall_asym}}
asymdhd <- function(){

}

difQpci<-function(x,y=NULL,q=seq(5,40,5)/100,xlab="Quantile",ylab="Group 1 minus Group 2",plotit=TRUE,alpha=.05,nboot=1000,SEED=TRUE,LINE=FALSE){
  #
  #  x can be a vector, in which case compare quantiels of distribution of data in x
  #  x can be a matrix with 2 columns, in which case analysis is done on dif=x[,1]=x[,2]
  #  y supplied, then do analysis of dif=x-y
  #
  #  Plot that provides perspective on the degree a distribution is symmetric about zero.
  #  This function plots the sum of q and 1-q quantiles. A 1-alpha confidence interval for the sum is indicated by a +
  #  If the distributions are symmetric
  #  the plot should be approximately a horizontal line. If in addition the median
  #  of the difference scores is zero, the horizontal line will intersect the y-axis at zero.
  #
  #  Similar to difQplot, only plots fewer quantiles by default and returns p-values for
  #  each quantile indicated by the argument q.
  #
  #  FWE is controlled via Hochberg's method, which was used to determine critical
  #  p-values based on the argument
  #  alpha.
  #
  #  Can alter the quantiles compared via the argument
  #  q
  #  q must be less than .5
  #
  #  LINE=TRUE. When plotting, a line connecting the estimates will be included.
  #
  x=as.matrix(x)
  if(is.null(y))dif=x
  if(ncol(x)>2)stop("Should be at most two groups")
  if(ncol(x)==2)dif=x[,1]-x[,2]
  if(!is.null(y))dif=x-y
  dif=elimna(dif)
  nv=length(dif)
  output=matrix(NA,ncol=8,nrow=length(q))
  dimnames(output)=list(NULL,c("quantile","Est_q","Est_1.minus.q","SUM","ci.low","ci.up","p_crit","p-value"))
  for(i in 1:length(q)){
    test=Dqdif(dif,q=q[i],plotit=FALSE,nboot=nboot,SEED=SEED)
    output[i,1]=q[i]
    output[i,2]=test$est.q
    output[i,3]=test$est.1.minus.q
    output[i,8]=test$p.value
    output[i,5]=test$conf.interval[1]
    output[i,6]=test$conf.interval[2]
  }
  temp=order(output[,8],decreasing=TRUE)
  zvec=alpha/c(1:length(q))
  output[temp,7]=zvec
  output <- data.frame(output)
  output$signif=rep("YES",nrow(output))
  for(i in 1:nrow(output)){
    if(output[temp[i],8]>output[temp[i],7])output$signif[temp[i]]="NO"
    if(output[temp[i],8]<=output[temp[i],7])break
  }
  output[,4]=output[,2]+output[,3]
  if(plotit){
    plot(rep(q,3),c(output[,4],output[,5],output[,6]),type="n",xlab=xlab,ylab=ylab)
    points(q,output[,6],pch="+")
    points(q,output[,5],pch="+")
    points(q,output[,4],pch="*")
    if(LINE)lines(q,output[,4],pch="*")
  }
  list(n=nv,output=output)
}

Dqdif<-function(x,y=NULL,q=.25,nboot=1000,plotit=TRUE,xlab="Group 1 - Group 2",SEED=TRUE,alpha=.05){
  #
  #  Compare two dependent groups by comparing the
  #  q and 1-q quantiles of the difference scores
  #
  # q should be < .5
  # if the groups do not differ, then the difference scores should be symmetric
  # about zero.
  # In particular, the sum of q and 1-q quantiles should be zero.
  #
  # q indicates the quantiles to be compared. By default, the .25 and .75 quantiles are used.
  #
  if(SEED)set.seed(2)
  if(q>=.5)stop("q should be less than .5")
  if(!is.null(y)){
    xy=elimna(cbind(x,y))
    dif=xy[,1]-xy[,2]
  }
  if(is.null(y))dif=elimna(x)
  n=length(dif)
  if(plotit)akerd(dif,xlab=xlab)
  bvec=NA
  data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
  for(ib in 1:nboot){
    bvec[ib]<-hd(dif[data[ib,]],q=q)+hd(dif[data[ib,]],q=1-q)
  }
  est1=hd(dif,q=q)
  est2=hd(dif,q=1-q)
  pv=mean(bvec<0)+.5*mean(bvec==0)
  p=2*min(c(pv,1-pv))
  low<-round((alpha/2)*nboot)+1
  up<-nboot-low
  sbvec=sort(bvec)
  ci=sbvec[low]
  ci[2]=sbvec[up]
  list(est.q=est1,est.1.minus.q=est2,conf.interval=ci,p.value=p)
}
