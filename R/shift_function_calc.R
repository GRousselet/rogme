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
#' @param data A data frame in tidy format. Column 1 describes the two groups;
#'   column 2 contains the values for each group. A properly formatted data
#'   frame can be created using \code{link{mkt2}}. Missing values are not
#'   allowed.
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
  #df <- na.omit(df) # remove NA
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
#' @param data A data frame in tidy format. Column 1 describes the two groups;
#'   column 2 contains the values for each group. A properly formatted data
#'   frame can be created using \code{link{mkt2}} or \code{link{mkt2d}}. Missing
#'   values are not allowed.
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
  #df <- na.omit(df) # remove NA
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
#' @param data A data frame in tidy format. Column 1 describes the two groups;
#'   column 2 contains the values for each group. A properly formatted data
#'   frame can be created using \code{link{mkt2}} or \code{link{mkt2d}}. Missing
#'   values are not allowed.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param q Quantiles to estimate - default = deciles 0.1:0.1:.9.
#' @param nboot Number of bootstrap samples - default = 2000
#' @param alpha Expected long-run type I error rate - default = 0.05
#' @param adj_ci Adjust confidence intervals & p values for multiple comparisons - default = TRUE
#' @return A data frame with one row per decile.
#' The columns are: \itemize{
#'   \item Column 1 = quantiles
#'   \item Column 2 = number of observations in group 1
#'   \item Column 3 = number of observations in group 2
#'   \item Column 4 = quantiles for group 1
#'   \item Column 5 = quantiles for group 2
#'   \item Column 6 = quantile differences (column 4 - column 5)
#'   \item Column 7 = lower bounds of the confidence intervals
#'   \item Column 8 = upper bounds of the confidence intervals
#'   \item Column 9 = critical p_values based on Hochberg's method
#'   \item Column 10 = p_values based on Hochberg's method
#'   \item Column 11 = significance 0/1
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
#' @export
shifthd_pbci <- function(data = df,
                         formula = obs ~ gr,
                         q = seq(.1, .9, .1),
                         nboot = 2000,
                         alpha = 0.05,
                         adj_ci = TRUE){
  # subset data
  out <- subset_data2(data, formula)
  #df <- na.omit(df) # remove NA
  x <- out$x
  y <- out$y
  nx <- length(x)
  ny <- length(y)
  gr_name1 <- out$gr_name1
  gr_name2 <- out$gr_name2
  # declare matrix of results
  output = matrix(0, nrow = length(q), ncol = 10)
  # confidence interval's boundaries
  low <- round((alpha / 2) * nboot) + 1
  up <- nboot - low
  # loop through quantiles
  for(i in 1:length(q)){
    output[i,1] = q[i]
    output[i,2] = nx
    output[i,3] = ny
    output[i,4] = hd(x, q = q[i])
    output[i,5] = hd(y, q = q[i])
    output[i,6] = output[i,4]-output[i,5]
    #  bootstrap
    datax <- matrix(sample(x, size = nx * nboot, replace = TRUE), nrow = nboot)
    datay <- matrix(sample(y, size = ny * nboot, replace = TRUE), nrow = nboot)
    bvecx <- apply(datax, 1, hd, q = q[i])
    bvecy <- apply(datay, 1, hd, q = q[i])
    bvec <- sort(bvecx - bvecy)
    temp <- sum(bvec < 0) / nboot + sum(bvec == 0) / (2 * nboot)
    output[i,7] = bvec[low] # ci_lower
    output[i,8] = bvec[up] # ci_upper
    output[i,10] = 2*(min(temp,1-temp)) # p_value
  }
  temp = order(output[,10], decreasing=TRUE)
  zvec = alpha / c(1:length(q))
  output[temp,9] = zvec # p_crit
  if(adj_ci){
    for(i in 1:length(q)){
      alpha = output[i,9]
      # confidence interval's boundaries
      low <- round((alpha / 2) * nboot) + 1
      up <- nboot - low
      #  bootstrap
      datax <- matrix(sample(x, size = nx * nboot, replace = TRUE), nrow = nboot)
      datay <- matrix(sample(y, size = ny * nboot, replace = TRUE), nrow = nboot)
      bvecx <- apply(datax, 1, hd, q = q[i])
      bvecy <- apply(datay, 1, hd, q = q[i])
      bvec <- sort(bvecx - bvecy)
      temp <- sum(bvec < 0) / nboot + sum(bvec == 0) / (2 * nboot)
      output[i,7] = bvec[low] # ci_lower
      output[i,8] = bvec[up] # ci_upper
      output[i,10] = 2*(min(temp,1-temp)) # p_value
    }
  }
  # make data frame
  out <- data.frame(output)
  names(out) <- c('q', 'n1', 'n2', gr_name1, gr_name2, 'difference',
                  'ci_lower', 'ci_upper', 'p_crit', 'p_value')
  # add sig column
  dplyr::mutate(out, sig = p_value <= p_crit)
  out
}

#' Shift function for two depend groups (pbci method)
#'
#' Compute a shift function for two dependent groups using the
#' Harrell-Davis quantile estimator in conjunction with a percentile
#' bootstrap approach.
#' Unlike \code{\link{shiftdhd}}: \itemize{
#' \item The confidence intervals are calculated using a percentile bootstrap of
#' the quantiles, instead of a percentile bootstrap of the standard error of
#' the difference of the quantiles.
#' \item The quantiles can be specified and are not limited to the deciles.
#' \item Tied values are allowed.
#' }
#' The confidence intervals are not corrected for multiple comparisons, the p
#' values are.
#' Examples of quantile sequences, from sparse to dense: \itemize{
#' \item \code{q = c(.25,.5,.75)}
#' \item \code{q = c(.1,.25,.5,.75,.9)}
#' \item \code{q = seq(.1, .9, .1)}
#' \item \code{q = seq(.05, .95, .05)}
#' }
#' @param data A data frame in tidy format. Column 1 describes the two groups;
#'   column 2 contains the values for each group. A properly formatted data
#'   frame can be created using \code{link{mkt2}} or \code{link{mkt2d}}. Missing
#'   values are not allowed.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param q Quantiles to estimate - default = deciles 0.1:0.1:.9.
#' @param nboot Number of bootstrap samples - default = 1000
#' @param alpha Expected long-run type I error rate - default = 0.05
#' @return A data frame with one row per decile.
#' The columns are: \itemize{
#'   \item Column 1 = quantiles
#'   \item Column 2 = number of observations
#'   \item Column 3 = quantiles for group 1
#'   \item Column 4 = quantiles for group 2
#'   \item Column 5 = quantile differences (column 3 - column 4)
#'   \item Column 6 = lower bounds of the confidence intervals
#'   \item Column 7 = upper bounds of the confidence intervals
#'   \item Column 8 = critical p_values based on Hochberg's method
#'   \item Column 9 = p_values
#'   \item Column 10 = significance 0/1
#'   }
#' @section Note:
#' Adaptation of Rand Wilcox's `Dqcomhd`, `bootdpci` & `rmmcppb` R functions
#' (\url{http://dornsife.usc.edu/labs/rwilcox/software/}).
#' From Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/}.
#'
#' @references
#' Wilcox, R.R. & Erceg-Hurn, D.M. (2012)
#' Comparing two dependent groups via quantiles.
#' J Appl Stat, 39, 2655-2664.
#' @seealso \code{hd}, \code{shiftdhd} for the pbse method for dependent
#'   groups, \code{shifthd_pbci} for dependent groups
#' @export
shiftdhd_pbci <- function(data = df,
                          formula = obs ~ gr,
                          q = seq(.1, .9, .1),
                          nboot = 1000,
                          alpha = 0.05){
  # subset data
  out <- subset_data2(data, formula)
  #df <- na.omit(df) # remove NA
  x <- out$x
  y <- out$y
  gr_name1 <- out$gr_name1
  gr_name2 <- out$gr_name2
  # declare matrix of results
  output = matrix(0, nrow = length(q), ncol = 9)
  # confidence interval's boundaries
  low <- round((alpha / 2) * nboot) + 1
  up <- nboot - low
  n <- length(x)
  # loop through quantiles
  for(i in 1:length(q)){
    output[i,1] = q[i]
    output[i,2] = n
    output[i,3] = hd(x, q = q[i])
    output[i,4] = hd(y, q = q[i])
    output[i,5] = output[i,3]-output[i,4]
    #  bootstrap
    bootsample <- matrix(sample(n, size = n * nboot, replace = TRUE), nrow = nboot)
    xmat <- matrix(x[bootsample], nrow = nboot, ncol = n)
    ymat <- matrix(y[bootsample], nrow = nboot, ncol = n)
    bvec <- apply(xmat, 1, hd, q = q[i]) - apply(ymat, 1, hd, q = q[i])
    bvec <- sort(bvec[,1] - bvec[,2])
    temp <- sum(bvec < 0) / nboot + sum(bvec == 0) / (2 * nboot)
    output[i,6] = bvec[low] # ci_lower
    output[i,7] = bvec[up] # ci_upper
    output[i,9] = 2*(min(temp,1-temp)) # p_value
  }
  temp = order(output[,9], decreasing=TRUE)
  zvec = alpha / c(1:length(q))
  output[temp,8] = zvec # p_crit
  # make data frame
  out <- data.frame(output)
  names(out) <- c('q', 'n', gr_name1, gr_name2, 'difference',
    'ci_lower', 'ci_upper', 'p_crit', 'p_value')
  # add sig column
  dplyr::mutate(out, sig = p_value <= p_crit)
  out
}
