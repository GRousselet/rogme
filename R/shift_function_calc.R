#' Shift function for two independent groups (pbse method)
#'
#' \code{shifthd} returns a shift function for two independent groups.
#' \code{shifthd} \strong{works exclusively with deciles and with alpha = 0.05},
#' so that only 95\% confidence intervals can be computed To use other quantiles
#' and other alpha levels, see \code{\link{shifthd_pbci}}.
#'
#' A shift function shows the difference between the quantiles of two groups as
#' a function of the quantiles of one group. For inferences, the function
#' returns a confidence interval for each quantile difference. The confidence
#' intervals are computed using a percentile bootstrap estimation of the
#' standard error of the difference between quantiles. The confidence intervals
#' are corrected for multiple comparisons so that the simultaneous probability
#' coverage is .95.
#'
#' The deciles are estimated using the Harrell-Davis quantile estimator of the
#' qth quantile is used. The default number of bootstrap samples is nboot = 200.
#' Independent bootstrap samples are used for each quantile and for each group.
#'
#' @section Note: Modified from Rallfun-v32.txt - see
#'   \url{https://github.com/nicebread/WRS/}.
#'
#' @references Harrell, F.E. & Davis, C.E. (1982) A new distribution-free
#' quantile estimator. Biometrika, 69, 635-640. Wilcox, R.R. (1995) Comparing
#' Two Independent Groups Via Multiple Quantiles. Journal of the Royal
#' Statistical Society. Series D (The Statistician), 44, 91-99.
#' @param df A data frame in tidy format. Column 1 describes the two groups;
#'   column 2 contains the values for each group. A properly formatted data
#'   frame can be created using \code{link{mkt2}}.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param nboot The number of bootstrap samples - default = 200.
#' @return A data frame with one row per decile. The columns are: \itemize{
#'   \item Column 1 = deciles for group 1 \item Column 2 = deciles for group 2
#'   \item Column 3 = differences (column 1 - column 2) \item Column 4 = lower
#'   bounds of the confidence intervals \item Column 5 = upper bounds of the
#'   confidence intervals }
#'
#' @seealso \code{\link{shiftdhd}} for dependent groups.
#'
#' @examples
#' out <- shifthd(df) # use the default parameters
#' out <- shifthd(df, nboot=500) # specify the number of bootstrap samples
#' out <- shifthd(df, ord(2,1)) # specify the order of the groups
#'
#' @export
shifthd <- function(data = df, formula = obs ~ gr, nboot = 200){
  # subset data
  out <- subset_data2(data, formula)
  df <- na.omit(df) # remove NA
  x <- out$x
  y <- out$y
  gr_name1 <- out$gr_name1
  gr_name2 <- out$gr_name2
  # factor to correct for multiple comparisons
  crit <- 80.1 / (min(length(x), length(y)))^2 + 2.73
  m <- matrix(0,9,5) # declare matrix of results
  # decile loop
  for (d in 1:9){
    q <- d/10
    # group 1
    data <- matrix(sample(x, size = length(x) * nboot, replace = TRUE),
                   nrow = nboot) # bootstrap samples
    bvec <- apply(data, 1, hd, q)
    se.x <- var(bvec)
    # group 2
    data <- matrix(sample(y, size = length(y) * nboot, replace = TRUE),
                   nrow = nboot) # bootstrap samples
    bvec <- apply(data, 1, hd, q)
    se.y <- var(bvec)
    m[d,1] <- hd(x,q)
    m[d,2] <- hd(y,q)
    m[d,3] <- m[d,1] - m[d,2]
    m[d,4] <- m[d,3] - crit * sqrt(se.x + se.y)
    m[d,5] <- m[d,3] + crit * sqrt(se.x + se.y)
  }
  out <- data.frame(m)
  names(out) <- c(gr_name1,gr_name2, 'difference', 'ci_lower', 'ci_upper')
  out
}


#' Shift function for two dependent groups (pbse methods)
#'
#' \code{shiftdhd} returns a shift function for two dependent groups.
#' \code{shiftdhd} \strong{works exclusively with deciles and with alpha =
#' 0.05}, so that only 95\% confidence intervals can be computed To use other
#' quantiles and other alpha levels, see \code{\link{shiftdhd_pbci}}.
#'
#' A shift function shows the difference between the quantiles of two groups as
#' a function of the quantiles of one group. For inferences, the function
#' returns a confidence interval for each quantile difference. The confidence
#' intervals are computed using a percentile bootstrap estimation of the
#' standard error of the difference between quantiles. The confidence intervals
#' are corrected for multiple comparisons so that the simultaneous probability
#' coverage is .95.
#'
#' The deciles are estimated using the Harrell-Davis quantile estimator of the
#' qth quantile is used. The default number of bootstrap samples is nboot = 200.
#' The same bootstrap samples are used for each quantile and for each group.
#'
#' @section Note: Modified from Rallfun-v32.txt - see
#'   \url{https://github.com/nicebread/WRS/}.
#'
#' @references Harrell, F.E. & Davis, C.E. (1982) A new distribution-free
#'   quantile estimator. Biometrika, 69, 635-640.
#'   Wilcox, R.R. (2012) Introduction to robust estimation and hypothesis
#'   testing. Academic Press, San Diego, CA.
#'
#' @param df A data frame in tidy format. Column 1 describes the two groups;
#'   column 2 contains the values for each group. A properly formatted data
#'   frame can be created using \code{link{mkt2}} or \code{link{mkt2d}}.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param nboot The number of bootstrap samples - default = 200.
#' @return A data frame with one row per decile. The columns are: \itemize{
#'   \item Column 1 = deciles for group 1 \item Column 2 = deciles for group 2
#'   \item Column 3 = differences (column 1 - column 2) \item Column 4 = lower
#'   bounds of the confidence intervals \item Column 5 = upper bounds of the
#'   confidence intervals }
#'
#' @seealso \code{\link{shifthd}} for independent groups.
#'
#' @examples
#' out <- shiftdhd(df, formula = obs ~ gr) # use the default parameters
#' out <- shiftdhd(df, formula = obs ~ gr, nboot = 500) # specify the number of bootstrap samples
#'
#' @export
shiftdhd <- function(data = df, formula = obs ~ gr, nboot = 200){
  # subset data
  out <- subset_data2(data, formula)
  df <- na.omit(df) # remove NA
  x <- out$x
  y <- out$y
  gr_name1 <- out$gr_name1
  gr_name2 <- out$gr_name2
  # factor to correct for multiple comparisons
  crit <- 37 / length(x)^(1.4) + 2.75
  m <- matrix(0,9,5) # declare matrix of results
  data <- matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
  xmat <- matrix(x[data],nrow=nboot,ncol=length(x))
  ymat <- matrix(y[data],nrow=nboot,ncol=length(x))
  for (d in 1:9){
    q <- d/10
    bvec <- apply(xmat, 1, hd, q) - apply(ymat, 1, hd, q)
    se <- sqrt(var(bvec))
    m[d,1]=hd(x,q)
    m[d,2]=hd(y,q)
    m[d,3]<-m[d,1]-m[d,2]
    m[d,4]<-m[d,3]-crit*se
    m[d,5]<-m[d,3]+crit*se
  }
  out <- data.frame(m)
  names(out) <- c(gr_name1, gr_name2, 'difference', 'ci_lower', 'ci_upper')
  out
}

# ===============================================================
#' Shift function for two independent groups (pbci method)
#'
#' Compute a shift function for two independent groups using the
#' Harrell-Davis quantile estimator in conjunction with a percentile
#' bootstrap approach.
#' Unlike \code{\link{shifthd}}: \itemize{
#' \item The confidence intervals are calculated using a percentile bootstrap of
#' the quantiles, instead of a percentile bootstrap of the standard error of
#' the difference of the quantiles.
#' \item The quantiles can be specified and are not limited to the deciles.
#' \item Tied values are allowed.
#' }
#' Examples of quantile sequences, from sparse to dense: \itemize{
#' \item \code{q = c(.25,.5,.75)}
#' \item \code{q = c(.1,.25,.5,.75,.9)}
#' \item \code{q = seq(.1, .9, .1)}
#' \item \code{q = seq(.05, .95, .05)}
#' }
#' @param x,y Vectors; no missing values are allowed
#' @param q Quantiles to estimate - default = deciles 0.1:0.1:.9.
#' Other suggestions
#' @param nboot Number of bootstrap samples - default = 2000
#' @param alpha Expected long-run type I error rate - default = 0.05
#' @return A data frame with one row per decile.
#' The columns are: \itemize{
#'   \item Column 1 = quantiles for group 1
#'   \item Column 2 = quantiles for group 2
#'   \item Column 3 = differences (column 1 - column 2)
#'   \item Column 4 = lower bounds of the confidence intervals
#'   \item Column 5 = upper bounds of the confidence intervals
#'   \item Column 6 = p-values based on Hochberg's method
#'   \item Column 7 = critical p-values based on Hochberg's method
#'   }
#' @section Note:
#' Adaptation of Rand Wilcox's `qcomhd` & `pb2gen` R functions
#' (\url{http://dornsife.usc.edu/labs/rwilcox/software/}).
#' From Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/}.
#'
#' @references
#' Wilcox, R.R., Erceg-Hurn, D.M., Clark, F. & Carlson, M. (2014)
#' Comparing two independent groups via the lower and upper quantiles.
#' J Stat Comput Sim, 84, 1543-1551.
#' @seealso \code{hd}, \code{shifthd} for the pbse method for independent
#'   groups, \code{shiftdhd_pbci} for dependent groups
shifthd_pbci <- function(data = df,
                         formula = obs ~ gr,
                         q = seq(.1, .9, .1),
                         nboot = 2000,
                         alpha = 0.05,
                         adj_ci = TRUE){
  # subset data
  out <- subset_data2(data, formula)
  df <- na.omit(df) # remove NA
  x <- out$x
  y <- out$y
  gr_name1 <- out$gr_name1
  gr_name2 <- out$gr_name2
  # declare matrix of results
  m = matrix(0, nrow = length(q), ncol = 9)
  # loop through quantiles
  for(i in 1:length(q)){
    output[i,1]=q[i]
    output[i,2]=length(elimna(x))
    output[i,3]=length(elimna(y))
    output[i,4]=hd(x,q=q[i])
    output[i,5]=hd(y,q=q[i])
    output[i,6]=output[i,4]-output[i,5]
    temp=pb2gen(x,y,nboot=nboot,est=hd,q=q[i],SEED=FALSE,alpha=alpha,pr=FALSE)
    output[i,7]=temp$ci[1]
    output[i,8]=temp$ci[2]
    output[i,10]=temp$p.value
  }
  temp=order(output[,10],decreasing=TRUE)
  zvec=alpha/c(1:length(q))
  output[temp,9]=zvec
  if(ADJ.CI){
    for(i in 1:length(q)){
      temp=pb2gen(x,y,nboot=nboot,est=hd,q=q[i],SEED=FALSE,alpha=output[i,9],pr=FALSE)
      output[i,7]=temp$ci[1]
      output[i,8]=temp$ci[2]
      output[i,10]=temp$p.value
    }
  }
  # make data frame
  out <- data.frame(m)
  names(out) <- c('q', gr_name1, gr_name2, 'difference',
                  'ci_lower', 'ci_upper', 'p_crit', 'p-value')
  # add sig column
  dplyr::mutate(out, sig = p-value <= p_crit)
  out
}

pb2gen<-function(x,y,alpha=.05,nboot=2000,est=onestep,SEED=TRUE,pr=FALSE,...){
  #
  #   Compute a bootstrap confidence interval for the
  #   the difference between any two parameters corresponding to
  #   independent groups.
  #   By default, M-estimators are compared.
  #   Setting est=mean, for example, will result in a percentile
  #   bootstrap confidence interval for the difference between means.
  #   Setting est=onestep will compare M-estimators of location.
  #   The default number of bootstrap samples is nboot=2000
  #
  x<-x[!is.na(x)] # Remove any missing values in x
  y<-y[!is.na(y)] # Remove any missing values in y
  if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
  bvecx<-apply(datax,1,est,...)
  bvecy<-apply(datay,1,est,...)
  bvec<-sort(bvecx-bvecy)
  low<-round((alpha/2)*nboot)+1
  up<-nboot-low
  temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
  sig.level<-2*(min(temp,1-temp))
  se<-var(bvec)
  list(est.1=est(x,...),est.2=est(y,...),est.dif=est(x,...)-est(y,...),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
}


#' Shift function for two depend groups (pbci methods)
#'
#' This version of the shift function
#' Computes a shift function for two dependent groups by comparing
#' the quantiles of the marginal distributions using the
#' Harrell-Davis quantile estimator in conjunction with a percentile
#' bootstrap approach.
#'
#' INPUTS:
#' - x & y are vectors of the same length without missing values
#' - q = quantiles to estimate - default = 0.1:0.1:.9 (deciles)
#'- nboot = number of bootstrap samples - default = 2000
#' - alpha = expected long-run type I error rate - default = 0.05
#' - plotit = 1 to get a figure; 0 otherwise by default
#'
#' OUTPUTS:
#'  - xd & yd = vectors of quantiles
#' - delta = vector of differences between quantiles (x-y)
#' - deltaCI = matrix quantiles x low/high bounds of the confidence
#' intervals of the quantile differences
#'
#' Unlike shiftdhd:
#' - the confidence intervals are not corrected for multiple comparisons
#' (the R version provides corrected critical p values - here i've decided
#' not to provide p values at all - the goal is to understand how
#' distributions differ, not to make binary decisions)
#' - the confidence intervals are calculated using a percentile bootstrap of
#'   the quantiles, instead of a percentile bootstrap of the standard error of the quantiles
#' - the quantiles to compare can be specified and are not limited to the
#'   deciles
#' - Tied values are allowed
#'
#' Unlike Rand Wilcox's Dqcomhd R function, no p value is returned.
#'  Extensive experience suggests humans cannot be trusted with p values.
#'
#' @section Note:
#' Adaptation of Rand Wilcox's Dqcomhd, bootdpci & rmmcppb R functions,
#' http://dornsife.usc.edu/labs/rwilcox/software/
#' From Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/}
#' @seealso HD, SHIFTDHD, SHIFTHD_PBCI
#' @references
#' Wilcox, R.R. & Erceg-Hurn, D.M. (2012)
#' Comparing two dependent groups via quantiles.
#' J Appl Stat, 39, 2655-2664.
shiftdhd_pbci <- function(){

}

# Dqcomhd<-function(x,y,q=c(1:9)/10,nboot=1000,plotit=TRUE,SEED=TRUE,xlab="Group 1",ylab="Est.1-Est.2",na.rm=TRUE,alpha=.05){
# #
# # Compare the quantiles of the marginal distributions associated with  two dependent groups
# # via hd estimator. Tied values are allowed.
# # When comparing lower or upper quartiles, both power and the probability of Type I error
# # compare well to other methods have been derived.
# #
# #  x: data for group 1
# #  y: data for group 2
# #  q: the quantiles to be compared
# #  nboot: Number of bootstrap samples
# #
# #
# if(SEED)set.seed(2)
#
# xy=elimna(cbind(x,y))
# x=xy[,1]
# y=xy[,2]
#
# pv=NULL
# output=matrix(NA,nrow=length(q),ncol=10)
# dimnames(output)<-list(NULL,c("q","n1","n2","est.1","est.2","est.1_minus_est.2","ci.low","ci.up","p_crit","p-value"))
# for(i in 1:length(q)){
# output[i,1]=q[i]
# output[i,2]=length(elimna(x))
# output[i,3]=length(elimna(y))
# output[i,4]=hd(x,q=q[i])
# output[i,5]=hd(y,q=q[i])
# output[i,6]=output[i,4]-output[i,5]
#
# temp=bootdpci(x,y,est=hd,q=q[i],dif=FALSE,plotit=FALSE,pr=FALSE,nboot=nboot,alpha=alpha,SEED=FALSE)
# output[i,7]=temp$output[1,5]
# output[i,8]=temp$output[1,6]
# output[i,10]=temp$output[1,3]
#
# }
# }
# if(plotit){
# xax=rep(output[,4],3)
# yax=c(output[,6],output[,7],output[,8])
# plot(xax,yax,xlab=xlab,ylab=ylab,type="n")
# points(output[,4],output[,6],pch="*")
# lines(output[,4],output[,6])
# points(output[,4],output[,7],pch="+")
# points(output[,4],output[,8],pch="+")
# }
# temp=order(output[,10],decreasing=TRUE)
# zvec=alpha/c(1:length(q))
# output[temp,9]=zvec
# output <- data.frame(output)
# output$signif=rep("YES",nrow(output))
# for(i in 1:nrow(output)){
# if(output[temp[i],10]>output[temp[i],9])output$signif[temp[i]]="NO"
# if(output[temp[i],10]<=output[temp[i],9])break
# }
# output
# }
