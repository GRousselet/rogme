#' Shift function for two independent groups (pbse method)
#'
#' \code{shifthd} returns a shift function for two independent groups or
#' multiple shift functions for pairs of independent groups. \code{shifthd}
#' \strong{works exclusively with deciles and with alpha = 0.05}, so that only
#' 95\% confidence intervals can be computed. To use other quantiles
#' and other alpha levels, see \code{\link{shifthd_pbci}}. Plot the shift
#' function using \code{plot_sf}.
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
#' qth quantile. The default number of bootstrap samples is nboot = 200.
#' Independent bootstrap samples are used for each quantile and for each group.
#'
#' @section Note: Modified from Rallfun-v32.txt - see
#'   \url{https://github.com/nicebread/WRS/}.
#'
#' @references Harrell, F.E. & Davis, C.E. (1982) A new distribution-free
#' quantile estimator. Biometrika, 69, 635-640.
#'
#' Wilcox, R.R. (1995) Comparing Two Independent Groups Via Multiple Quantiles.
#' Journal of the Royal Statistical Society. Series D (The Statistician), 44,
#' 91-99.
#'
#' @param data A data frame in long format. One column is a factor describing the groups;
#'   another column contains the values/observations for each group. A properly formatted data
#'   frame can be created using \code{\link{mkt2}}. Missing values are not
#'   allowed.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param nboot The number of bootstrap samples - default = 200.
#' @param todo A list of comparisons to perform - default = NULL.
#' @param doall Set to TRUE to compute all comparisons - default = FALSE. Not
#'   executed if a \code{todo} list is provided.
#'
#' @return A list of data frames, one data frame per comparison. Each data frame
#'   has one row per decile. The columns are: \itemize{ \item Column 1 = deciles
#'   for group 1 \item Column 2 = deciles for group 2
#'   \item Column 3 = differences (column 1 - column 2) \item Column 4 = lower
#'   bounds of the confidence intervals \item Column 5 = upper bounds of the
#'   confidence intervals }
#'
#' @seealso \code{\link{shiftdhd}} for dependent groups.
#'
#' \code{\link{plot_sf}} to plot the results.
#'
#' @examples
#' set.seed(21) # generate data
#' n <- 100 # sample size
#' df <- tibble(gr = c(rep("group1",n),rep("group2",n),rep("group3",n)),
#'              obs= c(rnorm(n)+6, rnorm(n)+4, rnorm(n)*1.5+6)) # make tibble
#'
#' out <- shifthd(df, obs ~ gr) # use the default parameters
#' out <- shifthd(df, obs ~ gr, nboot = 500) # specify the number of bootstrap samples
#' out <- shifthd(df, obs ~ gr, todo = list(c("g1","g2"),c("g3","g1"))) # specify list of comparisons
#' out <- shifthd(df, doall = TRUE) # compute all comparisons
#'
#' @export
shifthd <- function(data = df,
                    formula = obs ~ gr,
                    nboot = 200,
                    todo = NULL,
                    doall = FALSE){
  # subset data
  subf <- subset_formula(data, formula)
  if (length(todo)==0) { # no comparison is specified
    if (doall == FALSE) { # do not perform all comparisons
      if (length(subf$gr_names) > 2) {
        warning(paste0("Parameter column ",subf$param_col_name," contains more than 2 levels. The shift function is computed based on the first 2 levels: ",subf$gr_names[1], " vs. ",subf$gr_names[2]))
      }
      todo <- list(subf$gr_names[1:2])
    }
    if (doall == TRUE) { # perform all comparisons
      todo <- lapply(apply(combn(subf$gr_names, 2),2,list),unlist)
    }
  }
  out <- vector("list", length(todo)) # declare list of shift functions
  for(comp in 1:length(todo)){ # for each comparison
    x <- data[data[[subf$param_col_name]] == todo[[comp]][1], subf$obs_col_name][[1]]
    y <- data[data[[subf$param_col_name]] == todo[[comp]][2], subf$obs_col_name][[1]]
    gr_name1 <- todo[[comp]][1]
    gr_name2 <- todo[[comp]][2]
    # factor to correct for multiple comparisons
    crit <- 80.1 / (min(length(x), length(y)))^2 + 2.73
    m <- matrix(0,9,5) # declare matrix of results
    # decile loop
    for (d in 1:9){
      q <- d/10
      # group 1
      bootsamp <- matrix(sample(x, size = length(x) * nboot, replace = TRUE),
        nrow = nboot) # bootstrap samples
      bvec <- apply(bootsamp, 1, hd, q)
      se.x <- var(bvec)
      # group 2
      bootsamp <- matrix(sample(y, size = length(y) * nboot, replace = TRUE),
        nrow = nboot) # bootstrap samples
      bvec <- apply(bootsamp, 1, hd, q)
      se.y <- var(bvec)
      m[d,1] <- hd(x,q)
      m[d,2] <- hd(y,q)
      m[d,3] <- m[d,1] - m[d,2]
      m[d,4] <- m[d,3] - crit * sqrt(se.x + se.y)
      m[d,5] <- m[d,3] + crit * sqrt(se.x + se.y)
    }
    tmp <- data.frame(m)
    names(tmp) <- c(gr_name1, gr_name2, 'difference', 'ci_lower', 'ci_upper')
    out[[comp]] <- tmp
    names(out)[comp] <- paste0(gr_name1, " - ",gr_name2)
  }
  out
}


#' Shift function for two dependent groups (pbse methods)
#'
#' \code{shiftdhd} returns a shift function for two dependent groups or
#' multiple shift functions for pairs of dependent groups.
#' \code{shiftdhd} \strong{works exclusively with deciles and with alpha =
#' 0.05}, so that only 95\% confidence intervals can be computed. To use other
#' quantiles and other alpha levels, see \code{\link{shiftdhd_pbci}}. Plot the shift
#' function using \code{plot_sf}.
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
#' qth quantile. The default number of bootstrap samples is nboot = 200.
#' The same bootstrap samples are used for each quantile and for each group.
#'
#' @section Note: Modified from Rallfun-v32.txt - see
#'   \url{https://github.com/nicebread/WRS/}.
#'
#' @references Harrell, F.E. & Davis, C.E. (1982) A new distribution-free
#'   quantile estimator. Biometrika, 69, 635-640.
#'
#'   Wilcox, R.R. (2012) Introduction to robust estimation and hypothesis
#'   testing. Academic Press, San Diego, CA.
#'
#' @param data A data frame in long format. One column is a factor describing the conditions;
#'   another column contains the values/observations for each condition. A properly formatted data
#'   frame can be created using \code{link{mkt2}} or \code{link{mkt2d}}. Missing
#'   values are not allowed.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param nboot The number of bootstrap samples - default = 200.
#' @param todo A list of comparisons to perform - default = NULL.
#' @param doall Set to TRUE to compute all comparisons - default = FALSE. Not
#'   executed if a \code{todo} list is provided.
#'
#' @return A list of data frames, one data frame per comparison. Each data frame has one row per decile. The columns are: \itemize{
#'   \item Column 1 = deciles for group 1
#'   \item Column 2 = deciles for group 2
#'   \item Column 3 = differences (column 1 - column 2)
#'   \item Column 4 = lower bounds of the confidence intervals
#'   \item Column 5 = upper bounds of the confidence intervals }
#'
#' @seealso \code{\link{shifthd}} for independent groups.
#'
#' \code{\link{plot_sf}} to plot the results.
#'
#' @examples
#' set.seed(21) # generate data
#' n <- 100 # sample size
#' C1 <- rnorm(100)
#' df2 <- tibble(condition = factor(c(rep("C1",100),rep("C2",100),rep("C3",100))),
#'               data = c(C1+6, C1+rnorm(100)+4, C1+rnorm(100))) # make tibble
#'
#' out <- shiftdhd(df, obs ~ cond) # use the default parameters
#' out <- shiftdhd(df, obs ~ cond, nboot = 500) # specify the number of bootstrap samples
#' out <- shiftdhd(df, obs ~ cond, todo = list(c("C1","C2"),c("C3","C1"))) # specify list of comparisons
#' out <- shiftdhd(df, doall = TRUE) # compute all comparisons
#'
#' @export
shiftdhd <- function(data = df,
                     formula = obs ~ cond,
                     nboot = 200,
                     todo = NULL,
                     doall = FALSE){
  # subset data
  subf <- subset_formula(data, formula)
  # check all conditions have the same length
  if (length(unique(tapply(data[[subf$obs_col_name]], data[[subf$param_col_name]], length))) > 1) {
    stop("All conditions must have the same length")
  }
  if (length(todo)==0) { # no comparison is specified
    if (doall == FALSE) { # do not perform all comparisons
      if (length(subf$gr_names) > 2) {
        warning(paste0("Parameter column ",subf$param_col_name," contains more than 2 levels. The shift function is computed based on the first 2 levels: ",subf$gr_names[1], " vs. ",subf$gr_names[2]))
      }
      todo <- list(subf$gr_names[1:2])
    }
    if (doall == TRUE) { # perform all comparisons
      todo <- lapply(apply(combn(subf$gr_names, 2),2,list),unlist)
    }
  }
  out <- vector("list", length(todo)) # declare list of shift functions
  for(comp in 1:length(todo)){ # for each comparison
    x <- data[data[[subf$param_col_name]] == todo[[comp]][1], subf$obs_col_name][[1]]
    y <- data[data[[subf$param_col_name]] == todo[[comp]][2], subf$obs_col_name][[1]]
    gr_name1 <- todo[[comp]][1]
    gr_name2 <- todo[[comp]][2]
    # factor to correct for multiple comparisons
    crit <- 37 / length(x)^(1.4) + 2.75
    m <- matrix(0,9,5) # declare matrix of results
    # same bootstrap samples for all deciles and conditions (resample pairs)
    bootsamp <- matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
    xmat <- matrix(x[bootsamp],nrow=nboot,ncol=length(x))
    ymat <- matrix(y[bootsamp],nrow=nboot,ncol=length(x))
    for (d in 1:9){ # decile loop
      q <- d/10
      bvec <- apply(xmat, 1, hd, q) - apply(ymat, 1, hd, q)
      se <- sqrt(var(bvec))
      m[d,1]=hd(x,q)
      m[d,2]=hd(y,q)
      m[d,3]<-m[d,1]-m[d,2]
      m[d,4]<-m[d,3]-crit*se
      m[d,5]<-m[d,3]+crit*se
    }
    tmp <- data.frame(m)
    names(tmp) <- c(gr_name1, gr_name2, 'difference', 'ci_lower', 'ci_upper')
    out[[comp]] <- tmp
    names(out)[comp] <- paste0(gr_name1, " - ",gr_name2)
  }
  out
}

# ===============================================================
#' Shift function for two independent groups (pbci method)
#'
#' \code{shifthd_pbci} returns a shift function for two independent groups or
#' multiple shift functions for pairs of independent groups. It uses the
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
#' @param data A data frame in long format. One column is a factor describing the groups;
#'   another column contains the values/observations for each group. A properly formatted data
#'   frame can be created using \code{\link{mkt2}}. Missing values are not
#'   allowed.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param q Quantiles to estimate - default = deciles 0.1:0.1:.9.
#' @param nboot Number of bootstrap samples - default = 2000
#' @param alpha Expected long-run type I error rate - default = 0.05
#' @param adj_ci Adjust confidence intervals & p values for multiple comparisons - default = TRUE
#' @param todo A list of comparisons to perform - default = NULL.
#' @param doall Set to TRUE to compute all comparisons - default = FALSE. Not
#'   executed if a \code{todo} list is provided.
#' @return A list of data frames, one data frame per comparison. Each data frame
#'   has one row per decile. The columns are: \itemize{
#'   \item Column 1 = quantiles
#'   \item Column 2 = quantiles for group 1
#'   \item Column 3 = quantiles for group 2
#'   \item Column 4 = quantile differences (column 4 - column 5)
#'   \item Column 5 = lower bounds of the confidence intervals
#'   \item Column 6 = upper bounds of the confidence intervals
#'   \item Column 7 = critical p_values based on Hochberg's method
#'   \item Column 8 = p_values (based on Hochberg's method if adj_ci = TRUE)
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
#'
#' @seealso \code{\link{hd}}
#'
#' \code{\link{shifthd}} for the pbse method for independent groups,
#'
#' \code{\link{shiftdhd_pbci}} for dependent groups
#'
#' @export
shifthd_pbci <- function(data = df,
                         formula = obs ~ gr,
                         q = seq(.1, .9, .1),
                         nboot = 2000,
                         alpha = 0.05,
                         adj_ci = TRUE,
                         todo = NULL,
                         doall = FALSE){
  # subset data
  subf <- subset_formula(data, formula)
  if (length(todo)==0) { # no comparison is specified
    if (doall == FALSE) { # do not perform all comparisons
      if (length(subf$gr_names) > 2) {
        warning(paste0("Parameter column ",subf$param_col_name," contains more than 2 levels. The shift function is computed based on the first 2 levels: ",subf$gr_names[1], " vs. ",subf$gr_names[2]))
      }
      todo <- list(subf$gr_names[1:2])
    }
    if (doall == TRUE) { # perform all comparisons
      todo <- lapply(apply(combn(subf$gr_names, 2),2,list),unlist)
    }
  }
  # confidence interval's boundaries
  low <- round((alpha / 2) * nboot) + 1
  up <- nboot - low
  out <- vector("list", length(todo)) # declare list of shift functions
  for(comp in 1:length(todo)){ # for each comparison
    x <- data[data[[subf$param_col_name]] == todo[[comp]][1], subf$obs_col_name][[1]]
    y <- data[data[[subf$param_col_name]] == todo[[comp]][2], subf$obs_col_name][[1]]
    nx <- length(x)
    ny <- length(y)
    gr_name1 <- todo[[comp]][1]
    gr_name2 <- todo[[comp]][2]
    output <- matrix(0,length(q),8) # declare matrix of results
    # loop through quantiles
    for(i in 1:length(q)){
      output[i,1] = q[i]
      output[i,2] = hd(x, q = q[i])
      output[i,3] = hd(y, q = q[i])
      output[i,4] = output[i,2]-output[i,3]
      #  bootstrap
      datax <- matrix(sample(x, size = nx * nboot, replace = TRUE), nrow = nboot)
      datay <- matrix(sample(y, size = ny * nboot, replace = TRUE), nrow = nboot)
      bvecx <- apply(datax, 1, hd, q = q[i])
      bvecy <- apply(datay, 1, hd, q = q[i])
      bvec <- sort(bvecx - bvecy)
      temp <- sum(bvec < 0) / nboot + sum(bvec == 0) / (2 * nboot)
      output[i,5] = bvec[low] # ci_lower
      output[i,6] = bvec[up] # ci_upper
      output[i,8] = 2*(min(temp,1-temp)) # p_value
    }
    temp = order(output[,8], decreasing=TRUE)
    zvec = alpha / c(1:length(q))
    output[temp,7] = zvec # p_crit
    if(adj_ci){
      for(i in 1:length(q)){
        alpha = output[i,7]
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
        output[i,5] = bvec[low] # ci_lower
        output[i,6] = bvec[up] # ci_upper
        output[i,8] = 2*(min(temp,1-temp)) # p_value
      }
    }
    # make data frame
    tmp <- data.frame(output)
    names(tmp) <- c('q', gr_name1, gr_name2, 'difference',
      'ci_lower', 'ci_upper', 'p_crit', 'p_value')
    # add sig column
    # dplyr::mutate(tmp, sig = p_value <= p_crit)
    out[[comp]] <- tmp
    names(out)[comp] <- paste0(gr_name1, " - ",gr_name2)
  }
  out
}

#' Shift function for two depend groups (pbci method)
#'
#' \code{shiftdhd_pbci} returns a shift function for two independent groups or
#' multiple shift functions for pairs of independent groups. It uses the
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
#' @param data A data frame in long format. One column is a factor describing the groups;
#'   another column contains the values/observations for each group. A properly formatted data
#'   frame can be created using \code{\link{mkt2}}. Missing values are not
#'   allowed.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param q Quantiles to estimate - default = deciles 0.1:0.1:.9.
#' @param nboot Number of bootstrap samples - default = 1000
#' @param alpha Expected long-run type I error rate - default = 0.05
#' @param todo A list of comparisons to perform - default = NULL.
#' @param doall Set to TRUE to compute all comparisons - default = FALSE. Not
#'   executed if a \code{todo} list is provided.
#' @return A list of data frames, one data frame per comparison. Each data frame
#'   has one row per decile.
#' The columns are: \itemize{
#'   \item Column 1 = quantiles
#'   \item Column 2 = quantiles for group 1
#'   \item Column 3 = quantiles for group 2
#'   \item Column 4 = quantile differences (column 3 - column 4)
#'   \item Column 5 = lower bounds of the confidence intervals
#'   \item Column 6 = upper bounds of the confidence intervals
#'   \item Column 7 = critical p_values based on Hochberg's method
#'   \item Column 8 = p_values
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
#'
#' @seealso \code{\link{hd}}, \code{\link{shiftdhd}} for the pbse method for dependent
#'   groups, \code{\link{shifthd_pbci}} for independent groups
#'
#' @export
shiftdhd_pbci <- function(data = df,
                          formula = obs ~ gr,
                          q = seq(.1, .9, .1),
                          nboot = 1000,
                          alpha = 0.05,
                          todo = NULL,
                          doall = FALSE){
  # subset data
  subf <- subset_formula(data, formula)
  if (length(todo)==0) { # no comparison is specified
    if (doall == FALSE) { # do not perform all comparisons
      if (length(subf$gr_names) > 2) {
        warning(paste0("Parameter column ",subf$param_col_name," contains more than 2 levels. The shift function is computed based on the first 2 levels: ",subf$gr_names[1], " vs. ",subf$gr_names[2]))
      }
      todo <- list(subf$gr_names[1:2])
    }
    if (doall == TRUE) { # perform all comparisons
      todo <- lapply(apply(combn(subf$gr_names, 2),2,list),unlist)
    }
  }
  # confidence interval's boundaries
  low <- round((alpha / 2) * nboot) + 1
  up <- nboot - low
  out <- vector("list", length(todo)) # declare list of shift functions
  for(comp in 1:length(todo)){ # for each comparison
    x <- data[data[[subf$param_col_name]] == todo[[comp]][1], subf$obs_col_name][[1]]
    y <- data[data[[subf$param_col_name]] == todo[[comp]][2], subf$obs_col_name][[1]]
    n <- length(x)
    gr_name1 <- todo[[comp]][1]
    gr_name2 <- todo[[comp]][2]
    output <- matrix(0,length(q),8) # declare matrix of results
    # loop through quantiles
    for(i in 1:length(q)){
      output[i,1] = q[i]
      output[i,2] = hd(x, q = q[i])
      output[i,3] = hd(y, q = q[i])
      output[i,4] = output[i,2]-output[i,3]
      #  bootstrap
      bootsample <- matrix(sample(n, size = n * nboot, replace = TRUE), nrow = nboot)
      xmat <- matrix(x[bootsample], nrow = nboot, ncol = n)
      ymat <- matrix(y[bootsample], nrow = nboot, ncol = n)
      bvec <- apply(xmat, 1, hd, q = q[i]) - apply(ymat, 1, hd, q = q[i])
      bvec <- sort(bvec)
      temp <- sum(bvec < 0) / nboot + sum(bvec == 0) / (2 * nboot)
      output[i,5] = bvec[low] # ci_lower
      output[i,6] = bvec[up] # ci_upper
      output[i,8] = 2*(min(temp,1-temp)) # p_value
    }
    temp = order(output[,8], decreasing=TRUE)
    zvec = alpha / c(1:length(q))
    output[temp,7] = zvec # p_crit
    # make data frame
    tmp <- data.frame(output)
    names(tmp) <- c('q', gr_name1, gr_name2, 'difference',
      'ci_lower', 'ci_upper', 'p_crit', 'p_value')
    # add sig column
    # dplyr::mutate(tmp, sig = p_value <= p_crit)
    out[[comp]] <- tmp
    names(out)[comp] <- paste0(gr_name1, " - ",gr_name2)
  }
  out
}
