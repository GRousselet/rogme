#' Difference asymmetry function for two independent groups
#'
#' The difference asymmetry function provides perspective on the degree a distribution
#' is symmetric about zero, by quantifying the sum of q and 1-q quantiles of the
#' distribution of all pairwise differences between two independent groups x and y.
#' If the distributions of x and y are identical, then the
#' distribution of all pairwise differences is symmetric about zero.
#' In particular, the sum of q and 1-q quantiles should be zero.
#' If the distribution is symmetric the function should be approximately a horizontal line.
#' If in addition the median of the difference scores is zero, the horizontal line will
#' intersect the y-axis at zero.
#' Confidence intervals and p values are returned for each quantile sum.
#' The FWE is controlled via Hochberg's method, which is used to determine critical
#' p values based on the argument alpha.
#' To plot the results use \code{\link{plot_diff_asym}}.
#'
#' @param x,y Two vectors. Missing values are automatically removed.
#' @param q Quantiles to estimate - default = 0.05:0.05:.4 - must be <.5.
#' @param nboot Number of bootstrap samples - default = 1000.
#' @param alpha Expected long-run type I error rate - default = 0.05.
#'
#' @section Note:
#' This function combines Rand Wilcox's qwmwhd & cbmhd R functions,
#' from Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/},
#' http://dornsife.usc.edu/labs/rwilcox/software/.
#' @references
#' Wilcox, R.R., Erceg-Hurn, D.M., Clark, F. & Carlson, M. (2014)
#' Comparing two independent groups via the lower and upper quantiles.
#' J Stat Comput Sim, 84, 1543-1551.
#' @seealso \code{\link{hd}}, \code{\link{allpdiff}}, \code{\link{asymdhd}}
#' @export
asymhd <- function(x, y,
                   qseq = seq(5,40,5)/100,
                   alpha = .05,
                   nboot = 1000){
  x = elimna(x)
  y = elimna(y)
  n1 = length(x)
  n2 = length(y)
  m <- outer(x,y,FUN="-") # all pairwise differences
  output = matrix(NA,ncol=8,nrow=length(qseq))
  dimnames(output) = list(NULL,c("quantile","Est_q","Est_1.minus.q","SUM","ci.low","ci.up","p_crit","p-value"))
  for(i in 1:length(qseq)){
    q = qseq[i]
    output[i,1] = q
    if(q>=.5)stop("q should be less than .5")
    if(q<=0)stop("q should be greater than 0")
    q2 = 1-q
    est1 = hd(m,q)
    est2 = hd(m,q2)
    output[i,2] = est1
    output[i,3] = est2
    output[i,4] = est1 + est2
    # bootstrap samples
    data1 <- matrix(sample(n1,size=n1*nboot,replace=TRUE),nrow=nboot)
    data2 <- matrix(sample(n2,size=n2*nboot,replace=TRUE),nrow=nboot)
    bvec = NA
    for(Bi in 1:nboot){
      mb = outer(x[data1[Bi,]],y[data2[Bi,]],"-")
      bvec[Bi] = hd(mb,q) + hd(mb,q2)
    }
    # p value
    p = mean(bvec>0)+.5*mean(bvec==0)
    p = 2*min(c(p,1-p))
    output[i,8] = p
    # confidence interval
    sbv = sort(bvec)
    ilow <- round((alpha/2) * nboot)
    ihi <-  nboot - ilow
    ilow <- ilow+1
    ci = sbv[ilow]
    ci[2] = sbv[ihi]
    output[i,5] = ci[1]
    output[i,6] = ci[2]
  }
  temp = order(output[,8],decreasing=TRUE)
  zvec = alpha/c(1:length(qseq))
  output[temp,7] = zvec
  output <- data.frame(output)
  output$signif = rep("YES",nrow(output))
  for(i in 1:nrow(output)){
    if(output[temp[i],8] > output[temp[i],7])output$signif[temp[i]]="NO"
    if(output[temp[i],8] <= output[temp[i],7])break
  }
  list(n=c(n1,n2),output=output)
}

# ==========================================================================
#' Difference asymmetry function for two dependent groups or a single distribution
#'
#' \code{asymdhd} computes a difference asymmetry function for two dependent
#' groups or a single distribution. The difference asymmetry function provides
#' perspective on the degree a distribution is symmetric about zero, by
#' quantifying the sum of q and 1-q quantiles. If the groups do not differ, then
#' the difference scores should be symmetric about zero. In particular, the sum
#' of q and 1-q quantiles should be zero. If the distribution is symmetric the
#' function should be approximately a horizontal line. If in addition the median
#' of the difference scores is zero, the horizontal line will intersect the
#' y-axis at zero. Confidence intervals and p values are returned for each
#' quantile sum. The FWE is controlled via Hochberg's method, which is used to
#' determine critical p values based on the argument alpha. To plot the results
#' use \code{\link{plot_diff_asym}}.
#'
#' @param x A vector of pairwise differences.
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
#' @seealso \code{\link{hd}}, \code{\link{asymhd}}
#' @export
asymdhd <- function(x,
                    qseq = seq(5,40,5)/100,
                    alpha = .05,
                    nboot = 1000){
  dif = as.matrix(x)
  nv = length(dif)
  output = matrix(NA,ncol=8,nrow=length(qseq))
  dimnames(output) = list(NULL,c("quantile","Est_q","Est_1.minus.q","SUM","ci.low","ci.up","p_crit","p-value"))
  for(i in 1:length(qseq)){
    q = qseq[i]
    output[i,1] = q
    bvec = NA
    output[i,2] = hd(dif,q=q)
    output[i,3] = hd(dif,q=1-q)
    # bootstrap samples
    data <- matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
    for(ib in 1:nboot){
      bvec[ib] <- hd(dif[data[ib,]],q=q) + hd(dif[data[ib,]],q=1-q)
    }
    # p value
    pv = mean(bvec<0) + .5*mean(bvec==0)
    p = 2*min(c(pv,1-pv))
    output[i,8] = p
    # confidence interval
    low <- round((alpha/2)*nboot)+1
    up <- nboot-low
    sbvec = sort(bvec)
    ci = sbvec[low]
    ci[2] = sbvec[up]
    output[i,5] = ci[1]
    output[i,6] = ci[2]
  }
  temp = order(output[,8],decreasing=TRUE)
  zvec = alpha/c(1:length(q))
  output[temp,7] = zvec
  output <- data.frame(output)
  output$signif = rep("YES",nrow(output))
  for(i in 1:nrow(output)){
    if(output[temp[i],8] > output[temp[i],7]) output$signif[temp[i]]="NO"
    if(output[temp[i],8] <= output[temp[i],7]) break
  }
  output[,4] = output[,2] + output[,3]
  list(n = nv, output = output)
}

