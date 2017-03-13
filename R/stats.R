#' Compute deciles
#'
#' Estimate the deciles for the data in vector x using the Harrell-Davis
#' estimate of the qth quantile.
deciles <- function(x){
  xs<-sort(x)
  n<-length(x)
  vecx<-seq(along=x)
  xq<-0
  for (i in 1:9){
    q<-i/10
    m1<-(n+1)*q
    m2<-(n+1)*(1-q)
    wx<-pbeta(vecx/n,m1,m2)-pbeta((vecx-1)/n,m1,m2)  # W sub i values
    xq[i]<-sum(wx*xs)
  }
  xq
}

#' Estimates quantiles .1:.4 & .6:.9
#'
#' Estimate the deciles for the data in vector x using the Harrell-Davis
#' estimate of the qth quantile. Modified from \code{deciles} to return only
q1469 <- function(x){
  xq <- deciles(x)[c(seq(1,4),seq(6,9))]
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
  print(output)
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

# ==========================================================================
# Functions for the Kolmogorov-Smirnov test

#' Compute a Kolmogorov-Smirnov test
#'
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
