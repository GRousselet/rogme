#'  Shift function for two independent groups (pbse method)
#'
#'  \code{shifthd} returns a shift function for two independent groups.
#'  \code{shifthd} \strong{works exclusively with deciles and with alpha = 0.05}, so
#'  that only 95\% confidence intervals can be computed To use other quantiles
#'  and other alpha levels, see \code{\link{shifthd_pbci}}.
#'
#'  A shift function shows the difference between the quantiles of two groups as
#'  a function of the quantiles of one group. For inferences, the function
#'  returns a confidence interval for each quantile difference. The confidence
#'  intervals are computed using a percentile bootstrap estimation of the
#'  standard error of the difference between quantiles. The confidence intervals
#'  are corrected for multiple comparisons so that the simultaneous probability
#'  coverage is .95.
#'
#'  The deciles are estimated using the Harrell-Davis quantile estimator - see
#'  reference: of the qth quantile is used. The default number of bootstrap
#'  samples is nboot = 200.
#'
#' @section Note:
#' Modified from Rallfun-v31.txt - see \url{https://github.com/nicebread/WRS/}.
#'
#' @references Harrell, F.E. & Davis, C.E. (1982) A new distribution-free
#' quantile estimator. Biometrika, 69, 635-640. Wilcox, R.R. (1995) Comparing
#' Two Independent Groups Via Multiple Quantiles. Journal of the Royal
#' Statistical Society. Series D (The Statistician), 44, 91-99.
#' @param df A data frame in tidy format. Column 1 describes the two groups;
#'   column 2 contains the values for each group. A properly formatted
#'   data frame can be created using \code{link{mkt2}}.
#' @param nboot The number of bootstrap samples - default = 200.
#' @param ord Group order - default is to compute the difference between group 1
#'   and group 2 (group 1 - group 2).
#' @return A data frame with one row per decile. The columns are:
#' \itemize{
#'   \item Column 1 = deciles for group 1
#'   \item Column 2 = deciles for group 2
#'   \item Column 3 = differences (column 1 - column 2)
#'   \item Column 4 = lower bounds of the confidence intervals
#'   \item Column 5 = upper bounds of the confidence intervals
#' }
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
  names(out) <- c(gr_name1,gr_name2,'difference','ci_lower','ci_upper')
  out
}
