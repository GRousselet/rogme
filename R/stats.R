#' Compute quantiles
#'
#' Estimate the quantiles for the data in vector x using the Harrell-Davis
#' estimate of the qth quantile.
#' @param x A numeric vector of observations.
#' @param q A numeric vector of quantiles - default to the deciles.
#' @return A vector of quantile estimates.
#' @examples
#' x <- rnorm(100) # create vector
#' hdseq(x) # estimates of the deciles
#' hdseq(x, qseq = c(0.25, 0.5, 0.75)) # estimates of the quartiles
#' hdseq(x, qseq = c(seq(1,4),seq(6,9))) # estimates of the deciles except the median
#'
#' @export
hdseq <- function(x, qseq = seq(0.1,0.9,0.1)){
  xs <- sort(x)
  n <- length(x)
  vecx <- seq(along=x)
  xq <- vector(mode = "numeric", length = length(qseq))
  for (i in 1:length(qseq)){
    q <- qseq[i]
    m1 <- (n+1)*q
    m2 <- (n+1)*(1-q)
    wx <- pbeta(vecx/n,m1,m2) - pbeta((vecx-1)/n,m1,m2)  # W sub i values
    xq[i] <- sum(wx*xs)
  }
  xq
}

#' Percentile bootstrap confidence intervals of quantiles
#'
#' Compute percentile bootstrap confidence intervals of quantiles,
#' estimated using the Harrell-Davis estimator.
#' Default to deciles.
#' @section Note:
#' Combine functions \code{deciles} & \code{qcipb} from
#' Rallfun-v31.txt - see \url{https://github.com/nicebread/WRS}.
#' @export
quantiles_pbci <- function(x, q = seq(1,9) / 10, nboot = 2000, alpha = 0.05){
  low <- round((alpha/2)*nboot)
  up <- nboot-low
  low <- low+1
  x <- x[!is.na(x)]
  nq <- length(q)
  output = matrix(NA,ncol=4,nrow=nq)
  dimnames(output) = list(NULL,c("quantile","est_q","ci.low","ci.up"))
  bdata <- matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot) # bootstrap samples
  for (qi in 1:nq){
    output[qi,1] = q[qi]
    output[qi,2] = hd(x,q=q[qi])
    bvec <- apply(bdata,1,hd,q=q[qi])
    bvec <- sort(bvec)
    output[qi,3] = bvec[low]
    output[qi,4] = bvec[up]
  }
  output <- data.frame(output)
  # print(output)
  output
}

#' All pairwise differences
#'
#' Calculate all pairwise differences between 2 vectors.
#' @section Note:
#' Combine functions \code{loc2dif} & \code{wmwloc} from
#' Rallfun-v31.txt - see \url{https://github.com/nicebread/WRS}.
#' @param x,y Two vectors of same or different lengths.
#' @param paired If paired = TRUE and na.rm = TRUE, pairs are removed.
#' @param na.rm Missing values are removed if na.rm = TRUE.
#' @export
allpdiff <- function(x, y, paired = FALSE, na.rm = TRUE){
  if(na.rm){
    if(paired){
      m=elimna(cbind(x,y))
      x=m[,1]
      y=m[,2]
    } else {
      x<-x[!is.na(x)]
      y<-y[!is.na(y)]
    }
  }
  out <- as.vector(outer(x,y,FUN="-"))
  out
}

#' Confidence interval of the median of all pairwise differences
#'
#' Compute a confidence interval for the quantile (default to median)
#' of the distribution of D=X-Y (all pairwise differences),
#' where X and Y are two independent random variables.
#' The Harrell-Davis estimator is used.
#' If the distribution of X and Y are identical, then in particular the
#' median of the distribution of D=X-Y is approximately zero.
#' @param x,y Two vectors of continuous variables.
#' @param alpha Alpha level - default = 0.05, which means that a 95% confidence
#'   interval is computed.
#' @param q Quantile between 0 and 1 - default - 0.5, which means that the median of the distribution is estimated.
#' @param nboot Number of bootstrap samples - default = 600.
#' @return A list with variables: \code{q} \code{estimate}, \code{ci}, \code{n}
#' @section Note: Adaptation of \code{cbmhd} from Rallfun-v33.txt. See
#'   \url{https://github.com/nicebread/WRS/}
#' @seealso \code{\link{hd}} \code{\link{hdpbci}} \code{\link{cid}} \code{\link{allpdiff}}
#' @export
allpdiff_hdpbci <- function(x,y,alpha=.05,q=.5,nboot=600){
  if(q>=1)stop("q should be less than 1")
  if(q<=0)stop("q should be greater than 0")
  x<-x[!is.na(x)]
  y<-y[!is.na(y)]
  n1 = length(x)
  n2 = length(y)
  m <- outer(x, y, FUN = "-")
  n = n1*n2
  estimate = hd(m,q)
  # Bootstrap samples
  data1 <- matrix(sample(n1, size = n1*nboot, replace = TRUE), nrow = nboot)
  data2 <- matrix(sample(n2, size = n2*nboot, replace = TRUE), nrow = nboot)
  bvec=NA
  for(i in 1:nboot){
    mb = outer(x[data1[i,]],y[data2[i,]],"-")
    bvec [i] = hd(mb,q)
  }
  sbv = sort(bvec)
  ilow <- round((alpha/2) * nboot)
  ihi <- nboot - ilow
  ilow <- ilow+1
  ci = sbv[ilow]
  ci[2] = sbv[ihi]
  list(q = q, estimate = estimate, ci = ci, n = n)
}


#' Compare two independent groups using percentile bootstrap
#'
#' Compute a bootstrap confidence interval for the
#' the difference between any two parameters corresponding to
#' independent groups.
#' By default, the Harrell-Davis estimates of the median are compared.
#' @param x,y Two vectors of continuous variables.
#' @param alpha Alpha level - default = 0.05, which means that a 95% confidence
#'   interval is computed.
#' @param nboot Number of bootstrap samples - default = 2000.
#' @return A list with variables: \code{est.1} \code{est.2}, \code{est.diff} \code{ci}, \code{p.value}, \code{sq.se}, \code{n1}, \code{n2}
#' @section Note: Adaptation of \code{pb2gen} from Rallfun-v33.txt. See
#'   \url{https://github.com/nicebread/WRS/}
#' @seealso \code{\link{hd}}
#' @export
pb2gen <- function(x, y, alpha = .05, nboot = 2000, est = hd,...){
  x<-x[!is.na(x)] # Remove any missing values in x
  y<-y[!is.na(y)] # Remove any missing values in y
# Bootstrap samples
  datax <- matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  datay <- matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
  bvecx <- apply(datax,1,est,...)
  bvecy <- apply(datay,1,est,...)
  # Bootstrap confidence interval
  bvec <- sort(bvecx - bvecy)
  low <- round((alpha/2)*nboot)+1
  up <- nboot-low
  # P value
  temp <- sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
  sig.level<-2*(min(temp,1-temp))
  # Bootstrap estimate of the standard error of the difference
  se <- var(bvec)
  # Output
  list(est.1=est(x,...),est.2=est(y,...),est.dif=est(x,...)-est(y,...),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
}


#' Harrell-Davis quantile estimator
#'
#' Compute the Harrell-Davis estimate of the qth quantile.
#'
#' \code{hd} is a weighted average of all the order statistics. See
#' \url{https://garstats.wordpress.com/2016/06/09/the-harrell-davis-quantile-estimator/}
#' for more details.
#' @param x A vector of continuous observations.
#' @param q A quantile between 0 and 1 - default is 0.5, the median.
#' @param na.rm Set to TRUE to remove NA, FALSE otherwise - default is TRUE.
#' @examples
#' hd(x) # use default estimation of the median (q = 0.5)
#' hd(x, q = 0.25) # estimate the first quartile
#' hd(x, q = 0.75) # estimate the third quartile
#' @references
#' Harrell, F.E. & Davis, C.E. (1982) A new distribution-free quantile estimator. Biometrika, 69, 635-640.
#' @section Note: code from Rallfun-v33.txt. See
#'   \url{https://github.com/nicebread/WRS/}
#' @export
hd <- function(x,q=.5,na.rm=TRUE){
  if(na.rm) x <- x[complete.cases(x)]
  n <- length(x)
  m1 <- (n + 1) * q
  m2 <- (n + 1) * (1 - q)
  vec <- seq(along = x)
  w <- pbeta(vec / n, m1, m2) - pbeta((vec - 1)/n, m1, m2)  # weights
  y <- sort(x)
  hd <- sum(w * y) # weighted sum
  hd
}

#' Bootstrap confidence interval of the qth quantile
#'
#' Compute a percentile bootstrap confidence interval for the qth quantile via
#' the Harrell--Davis estimator. Appears to be best method when there are tied
#' values.
#' @param x A vector of continuous observations.
#' @param q A quantile between 0 and 1 - default = 0.5, meaning a confidence
#'   interval for the median is computed.
#' @param alpha Alpha level - default = 0.05, such that a 95% confidence
#'   interval is computed.
#' @param nboot Number of bootstrap samples - default = 2000.
#' @param nv Null value when computing a p-value.
#' @return A list with variables: \code{q} \code{estimate}, \code{ci}, \code{n},
#'   \code{p.value}
#' @section Note: Adaptation of \code{qcipb} from Rallfun-v33.txt. See
#'   \url{https://github.com/nicebread/WRS/}
#' @seealso \code{\link{hd}}
#' @export
hdpbci <- function(x, q = .5, alpha = .05, nboot = 2000, nv = 0){
  x = elimna(x)
  # Compute quantile
  estimate = hd(x,q=q)
  # Bootstrap estimates
  data <- matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  bvec <- apply(data,1,hd,q=q)
  # Bootstrap confidence intervals
  bvec <- sort(bvec)
  low <- round((alpha/2)*nboot)
  up <- nboot-low
  low <- low+1
  # P value
  pv = mean(bvec>nv) + .5*mean(bvec==nv)
  pv = 2*min(c(pv,1-pv))
  # Outputs
  list(q=q, estimate=estimate, ci=c(bvec[low],bvec[up]), n=length(x), p.value=pv)
}

# ==========================================================================
# Functions for the Kolmogorov-Smirnov test

#' Compute the Kolmogorov-Smirnov test statistic
#'
#' \code{ks} returns a list containing the value of the test statistic, the
#' approximate .05 critical value, and the exact significance level if sig=T.
#'
#' @param x,y Two vectors. Missing values are automatically removed.
#' @param w TRUE or FALSE. Set w to FALSE to use the standard version. Set w to
#'   TRUE to use a weighted version, which is more sensitive to differences
#'   occuring in the tails of the distributions.
#' @param sig TRUE to use the exact significance level. If there are ties, the
#'   reported significance level is exact when using the unweighted test, but
#'   for the weighted test the reported level is too high.
#' @seealso This function uses the functions ecdf, kstiesig, kssig and kswsig
#' @section Note:
#' From Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/}
#' @export
ks <- function(x, y, w = FALSE, sig = TRUE, alpha = .05){
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  w <- as.logical(w)
  sig <- as.logical(sig)
  tie <- logical(1)
  siglevel <- NA
  z <- sort(c(x,y))  # Pool and sort the observations
  tie = duplicated(z)
  v <- 1   # Initializes v
  for (i in 1:length(z))v[i]<-abs(ecdf(x,z[i])-ecdf(y,z[i]))
  ks<-max(v)
  #
  #crit<-1.36*sqrt((length(x)+length(y))/(length(x)*length(y))) # Approximate
  #                                                       .05 critical value
  crit=sqrt(0-log(alpha/2)*(length(x)+length(y))/(2*length(x)*length(y)))
  if(!w && sig && !tie)siglevel<-kssig(length(x),length(y),ks)
  if(!w && sig && tie)siglevel<-kstiesig(x,y,ks)
  if(w){
    crit<-(max(length(x),length(y))-5)*.48/95+2.58+abs(length(x)-length(y))*.44/95
    if(length(x)>100 || length(y)>100)warning(paste("When either sample size is
      greater than 100, the approximate critical value can be inaccurate. It is
      recommended that the exact significance level be computed."))
    for (i in 1:length(z)){
      temp<-(length(x)*ecdf(x,z[i])+length(y)*ecdf(y,z[i]))/length(z)
      temp<-temp*(1.-temp)
      v[i]<-v[i]/sqrt(temp)
    }
    v<-v[!is.na(v)]
    ks<-max(v)*sqrt(length(x)*length(y)/length(z))
    if(sig)siglevel<-kswsig(length(x),length(y),ks)
    if(tie && sig)
      warning(paste("Ties were detected. The reported significance level of the
        weighted Kolmogorov-Smirnov test statistic is not exact."))
  }
  list(test = ks, critval = crit, p.value = siglevel)
}

#'  Compute empirical cdf
#'  Compute the empirical cdf for data in x evaluated at val.
#'  That is, estimate P(X <= val)
ecdf <- function(x, val){
  ecdf <- length(x[x <= val]) / length(x)
  ecdf
}


#' Compute significance level of the Kolmogorov-Smirnov test statistic
#'
#' Compute significance level of the Kolmogorov-Smirnov test statistic
#' for the data in vectors x and y at val, the observed value of the
#' test statistic.
#' \code{kstiesig} allows ties among values.
kstiesig <- function(x, y, val){
  m<-length(x)
  n<-length(y)
  z<-c(x,y)
  z<-sort(z)
  cmat<-matrix(0,m+1,n+1)
  umat<-matrix(0,m+1,n+1)
  for (i in 0:m){
    for (j in 0:n){
      if(abs(i/m-j/n)<=val)cmat[i+1,j+1]<-1e0
      k<-i+j
      if(k > 0 && k<length(z) && z[k]==z[k+1])cmat[i+1,j+1]<-1
    }
  }
  for (i in 0:m){
    for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
    else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
  }
  term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
  kstiesig<-1.-umat[m+1,n+1]/exp(term)
  kstiesig
}

#' Compute significance level of the  Kolmogorov-Smirnov test statistic
#' m=sample size of first group
#' n=sample size of second group
#' val=observed value of test statistic
kssig<-function(m,n,val){
  cmat<-matrix(0,m+1,n+1)
  umat<-matrix(0,m+1,n+1)
  for (i in 0:m){
    for (j in 0:n)cmat[i+1,j+1]<-abs(i/m-j/n)
  }
  cmat<-ifelse(cmat<=val,1e0,0e0)
  for (i in 0:m){
    for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
    else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
  }
  term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
  kssig<-1.-umat[m+1,n+1]/exp(term)
  kssig
}

#' Compute significance level of the weighted
#' Kolmogorov-Smirnov test statistic
#'
#' @param m Sample size of first group
#' @param n Sample size of second group
#' @param val Observed value of test statistic
kswsig <- function(m, n, val){
  mpn<-m+n
  cmat<-matrix(0,m+1,n+1)
  umat<-matrix(0,m+1,n+1)
  for (i in 1:m-1){
    for (j in 1:n-1)cmat[i+1,j+1]<-abs(i/m-j/n)*sqrt(m*n/((i+j)*(1-(i+j)/mpn)))
  }
  cmat<-ifelse(cmat<=val,1,0)
  for (i in 0:m){
    for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
    else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
  }
  term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
  kswsig<-1.-umat[m+1,n+1]/exp(term)
  kswsig
}
# ==========================================================================

# ==========================================================================
# FUNCTIONS FOR CLIFF'S DELTA TEST

#' P value of Cliff's delta test
#'
#' Compute Cliff's analog of WMW test & associated p value
#'
#' \code{cidv2} returns a list containing the value of the test statistic, its
#' confidence interval, and its p value.
#'
#' @param x,y Two vectors. Missing values are automatically removed.
#' @param alpha Alpha significance level.
#' @seealso This function uses the function cid.
#' @section Note:
#' From Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/}
#' @references
#' Cliff, N. (1996) Ordinal methods for behavioral data analysis. Erlbaum, Mahwah, N.J.
#' Wilcox, R.R. (2011) Inferences about a Probabilistic Measure of Effect Size When Dealing with More Than Two Groups. Journal of Data Science, 9, 471-486.
#' @export
cidv2 <- function(x,y,alpha=.05){
  nullval <- 0
  ci <- cid(x,y,alpha=alpha)
  alph<-c(1:99)/100
  for(i in 1:99){
    irem<-i
    chkit <- cid(x,y,alpha=alph[i])
    if(chkit[[3]]>nullval || chkit[[4]]<nullval)break
  }
  p.value<-irem/100
  if(p.value<=.1){
    iup<-(irem+1)/100
    alph<-seq(.001,iup,.001)
    for(i in 1:length(alph)){
      p.value<-alph[i]
      chkit<-cid(x,y,alpha=alph[i])
      if(chkit[[3]]>nullval || chkit[[4]]<nullval)break
    }}
  if(p.value<=.001){
    alph<-seq(.0001,.001,.0001)
    for(i in 1:length(alph)){
      p.value<-alph[i]
      chkit<-cid(x,y,alpha=alph[i])
      if(chkit[[3]]>nullval || chkit[[4]]<nullval)break
    }}
  phat<-(1-ci$d)/2
  pci=c((1-ci$cu)/2,(1-ci$cl)/2)
  d.ci=c(ci$cl,ci$cu)
  dval=cid(x,y)$summary.dvals
  list(n1=length(elimna(x)),n2=length(elimna(y)),d.hat=ci$d,d.ci=d.ci,p.value=p.value,p.hat=phat,p.ci=pci,summary.dvals=dval)
}

#' Cliff's delta test
#'
#' Compute Cliff's analog of WMW test. Compute a confidence interval for delta
#' using the method in Cliff, 1996, p. 140, eq 5.12.
#' The null hypothesis is that for two independent group, P(X<Y)=P(X>Y).
#' This function reports a 1-alpha confidence interval for
#' P(X>Y)-P(X<Y).
#'
#' \code{cid} returns a list containing the value of the test statistic, its
#' confidence interval, and its p value.
#'
#' @param x,y Two vectors. Missing values are automatically removed.
#' @param alpha Alpha significance level.
#' @seealso To compute a p value, use \code{\link{cidv2}} \code{\link{allpdiff_hdpbci}}
#' @section Note:
#' From Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/}
#' @references
#' Cliff, N. (1996) Ordinal methods for behavioral data analysis. Erlbaum, Mahwah, N.J.
#' @export
cid <- function(x, y, alpha=.05){
  x<-x[!is.na(x)]
  y<-y[!is.na(y)]
  m<-outer(x,y,FUN="-")
  msave<-m
  m<-sign(m)
  d<-mean(m)
  phat<-(1-d)/2
  flag=T
  if(phat==0 || phat==1)flag=F
  q0<-sum(msave==0)/length(msave)
  qxly<-sum(msave<0)/length(msave)
  qxgy<-sum(msave>0)/length(msave)
  c.sum<-matrix(c(qxly,q0,qxgy),nrow=1,ncol=3)
  dimnames(c.sum)<-list(NULL,c("P(X<Y)","P(X=Y)","P(X>Y)"))
  if(flag){
    sigdih<-sum((m-d)^2)/(length(x)*length(y)-1)
    di<-NA
    for (i in 1:length(x))di[i]<-sum(x[i]>y)/length(y)-sum(x[i]<y)/length(y)
    dh<-NA
    for (i in 1:length(y))dh[i]<-sum(y[i]>x)/length(x)-sum(y[i]<x)/length(x)
    sdi<-var(di)
    sdh<-var(dh)
    sh<-((length(y)-1)*sdi+(length(x)-1)*sdh+sigdih)/(length(x)*length(y))
    zv<-qnorm(alpha/2)
    cu<-(d-d^3-zv*sqrt(sh)*sqrt((1-d^2)^2+zv^2*sh))/(1-d^2+zv^2*sh)
    cl<-(d-d^3+zv*sqrt(sh)*sqrt((1-d^2)^2+zv^2*sh))/(1-d^2+zv^2*sh)
  }
  if(!flag){
    sh=NULL
    nm=max(c(length(x),length(y)))
    if(phat==1)bci=binomci(nm,nm,alpha=alpha)
    if(phat==0)bci=binomci(0,nm,alpha=alpha)
  }
  if(flag)pci=c((1-cu)/2,(1-cl)/2)
  if(!flag){
    pci=bci$ci
    cl=1-2*pci[2]
    cu=1-2*pci[1]
  }
  list(n1=length(x),n2=length(y),cl=cl,cu=cu,d=d,sqse.d=sh,phat=phat,summary.dvals=c.sum,ci.p=pci)
}

binomci <- function(x=sum(y),nn=length(y),y=NULL,n=NA,alpha=.05){
  #  Compute a 1-alpha confidence interval for p, the probability of
  #  success for a binomial distribution, using Pratt's method
  #
  #  y is a vector of 1s and 0s.
  #  x is the number of successes observed among n trials
  #
  # From Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/}
  if(!is.null(y)){
    y=elimna(y)
    nn=length(y)
  }
  if(nn==1)stop("Something is wrong: number of observations is only 1")
  n<-nn
  if(x!=n && x!=0){
    z<-qnorm(1-alpha/2)
    A<-((x+1)/(n-x))^2
    B<-81*(x+1)*(n-x)-9*n-8
    C<-(0-3)*z*sqrt(9*(x+1)*(n-x)*(9*n+5-z^2)+n+1)
    D<-81*(x+1)^2-9*(x+1)*(2+z^2)+1
    E<-1+A*((B+C)/D)^3
    upper<-1/E
    A<-(x/(n-x-1))^2
    B<-81*x*(n-x-1)-9*n-8
    C<-3*z*sqrt(9*x*(n-x-1)*(9*n+5-z^2)+n+1)
    D<-81*x^2-9*x*(2+z^2)+1
    E<-1+A*((B+C)/D)^3
    lower<-1/E
  }
  if(x==0){
    lower<-0
    upper<-1-alpha^(1/n)
  }
  if(x==1){
    upper<-1-(alpha/2)^(1/n)
    lower<-1-(1-alpha/2)^(1/n)
  }
  if(x==n-1){
    lower<-(alpha/2)^(1/n)
    upper<-(1-alpha/2)^(1/n)
  }
  if(x==n){
    lower<-alpha^(1/n)
    upper<-1
  }
  phat<-x/n
  list(phat=phat,ci=c(lower,upper),n=n)
}
# ==========================================================================

#' Proportion of observations greater than a specified value
#'
#' For a vector x, return the proportion of observations greater than a specified value a: P(X>a).
#' Also return the complementary proportions P(X<a) and P(X=a).
#'
#' @param x A numeric vector
#' @param a Value for comparison  - default = 0
#' @return A list of 3 elements: \itemize{
#' \item \code{P(X>a)}
#' \item \code{P(X=a)}
#' \item \code{P(X<a)}
#' }
#'
#' @seealso adapted from \code{\link{cid}}
#' @section Note:
#' From Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/}
#' @export
pxgta <- function(x, a = 0){
  x<-x[!is.na(x)]
  if(!is.numeric(x)){
    stop("input x must be numeric")
  }
  if(!is.numeric(a)){
    stop("input a must be numeric")
  }
  pxlta <- sum(x < a) / length(x)
  pa <- sum(x==a) / length(x)
  pxgta <- sum(x > a) / length(x)
  out <- list(pxlta, pa, pxgta)
  names(out) <- c("P(X<a)","P(X=a)","P(X>a)")
  out
}

#' Proportion of observations in x greater than observations in y
#'
#' For two vectors x and y, return the proportion of observations in x greater
#' than observations in y: P(X>Y). Also return the complementary proportions
#' P(X<Y) and P(X=Y). The proportions are determined by computing all pariwise
#' differences between vectors.
#'
#' @param x A numeric vector
#' @param y A numeric vector
#' @return A list of 3 elements: \itemize{
#' \item \code{P(X>Y)}
#' \item \code{P(X=Y)}
#' \item \code{P(X<Y)}
#' }
#'
#' @seealso adapted from \code{\link{cid}}
#' @section Note:
#' From Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/}
#' @references
#' Cliff, N. (1996) Ordinal methods for behavioral data analysis. Erlbaum, Mahwah, N.J.
#' @export
pxgty <- function(x, y){
  x<-x[!is.na(x)]
  if(!is.numeric(x)){
    stop("input x must be numeric")
  }
  if(!is.numeric(y)){
    stop("input y must be numeric")
  }
  m <- outer(x, y, FUN="-")
  pxlty <- sum(m < 0) / length(m)
  p0 <- sum(m==0) / length(m)
  pxgty <- sum(m > 0) / length(m)
  out <- list(pxlty, p0, pxgty)
  names(out) <- c("P(X<Y)","P(X=0)","P(X>Y)")
  out
}
