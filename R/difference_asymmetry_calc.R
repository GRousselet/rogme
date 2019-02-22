#' Difference asymmetry function for two independent groups
#'
#' The difference asymmetry function provides perspective on the degree a
#' distribution is symmetric about zero, by quantifying the sum of q and 1-q
#' quantiles of the distribution of all pairwise differences between two
#' independent groups x and y. If the distributions of x and y are identical,
#' then the distribution of all pairwise differences is symmetric about zero. In
#' particular, the sum of q and 1-q quantiles should be zero. If the
#' distribution is symmetric the function should be approximately a horizontal
#' line. If in addition the median of the difference scores is zero, the
#' horizontal line will intersect the y-axis at zero. Confidence intervals and p
#' values are returned for each quantile sum. The FWE is controlled via
#' Hochberg's method, which is used to determine critical p values based on the
#' argument alpha. To plot the results use \code{\link{plot_diff_asym}}.
#'
#' @param data A data frame in long format. One column is a factor describing the groups;
#'   another column contains the values/observations for each group. A properly formatted data
#'   frame can be created using \code{\link{mkt2}}. Missing values are not
#'   allowed.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param qseq Quantiles to estimate - default = 0.05:0.05:.4 - must be <.5.
#' @param nboot Number of bootstrap samples - default = 1000.
#' @param alpha Expected long-run type I error rate - default = 0.05.
#' @param todo A list of comparisons to perform - default = NULL.
#' @param doall Set to TRUE to compute all comparisons - default = FALSE. Not
#'   executed if a \code{todo} list is provided.
#' @return A list of data frames, one data frame per comparison. Each data frame
#'   has one row per decile. The columns are: \itemize{
#'   \item Column 1 = quantiles
#'   \item Column 2 = quantiles of differences
#'   \item Column 3 = 1 - quantiles of differences
#'   \item Column 4 = sum of quantiles
#'   \item Column 5 = lower bounds of the confidence intervals
#'   \item Column 6 = upper bounds of the confidence intervals
#'   \item Column 7 = critical p_values based on Hochberg's method
#'   \item Column 8 = p_values
#'   }
#'
#' @section Note:
#' This function combines Rand Wilcox's qwmwhd & cbmhd R functions,
#' from Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/},
#' http://dornsife.usc.edu/labs/rwilcox/software/.
#'
#' @references
#' Wilcox, R.R., Erceg-Hurn, D.M., Clark, F. & Carlson, M. (2014)
#' Comparing two independent groups via the lower and upper quantiles.
#' J Stat Comput Sim, 84, 1543-1551.
#'
#' @seealso \code{\link{hd}} Harrell-Davis quantile estimator
#' \code{\link{allpdiff}} to compute all pairwise differences
#' \code{\link{asymdhd}} for dependent groups
#' \code{\link{plot_diff_asym}} to plot results
#'
#' @examples
#' set.seed(21) # generate data
#' n <- 100 # sample size
#' df <- tibble(gr = factor(c(rep("group1",n),rep("group2",n),rep("group3",n))),
#'              obs= c(rnorm(n)+6, rnorm(n)+4, rnorm(n)*1.5+6)) # make tibble
#'
#' out <- asymhd(df, obs ~ gr) # use the default parameters
#' out <- asymhd(df, obs ~ gr, alpha = .90) # specify alpha
#' out <- asymhd(df, obs ~ gr, nboot = 500) # specify the number of bootstrap samples
#' out <- asymhd(df, obs ~ gr, todo = list(c("group1","group2"),c("group3","group1"))) # specify list of comparisons
#' out <- asymhd(df, obs ~ gr, q = seq(.1,.4,.1)) # specify the quantiles
#' out <- asymhd(df, obs ~ gr, doall = TRUE) # compute all comparisons
#'
#' @export
asymhd <- function(data = df,
                   formula = obs ~ gr,
                   qseq = seq(5,40,5)/100,
                   alpha = .05,
                   nboot = 1000,
                   todo = NULL,
                   doall = FALSE){
  # check input is a data frame
  if(!is.data.frame(data)){
    stop("data must be a data frame")
  }
  # subset data
  subf <- subset_formula(data, formula)
  if (length(todo)==0) { # no comparison is specified
    if (doall == FALSE) { # do not perform all comparisons
      if (length(subf$gr_names) > 2) {
        warning(paste0("Parameter column ",subf$param_col_name," contains more than 2 levels. The difference asymmetry function is computed based on the first 2 levels: ",subf$gr_names[1], " vs. ",subf$gr_names[2]))
      }
      todo <- list(subf$gr_names[1:2])
    }
    if (doall == TRUE) { # perform all comparisons
      todo <- lapply(apply(combn(subf$gr_names, 2),2,list),unlist)
    }
  }
  # confidence interval's boundaries
  ilow <- round((alpha/2) * nboot)
  ihi <-  nboot - ilow
  ilow <- ilow+1
  out <- vector("list", length(todo)) # declare list of difference asymmetry functions
  for(comp in 1:length(todo)){ # for each comparison
    x <- data[data[[subf$param_col_name]] == todo[[comp]][1], subf$obs_col_name][[1]]
    y <- data[data[[subf$param_col_name]] == todo[[comp]][2], subf$obs_col_name][[1]]
    nx <- length(x)
    ny <- length(y)
    gr_name1 <- todo[[comp]][1]
    gr_name2 <- todo[[comp]][2]
    output <- matrix(0,nrow=length(qseq),ncol=8) # declare matrix of results
    m <- outer(x,y,FUN="-") # all pairwise differences
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
    bootdatax <- matrix(sample(nx,size=nx*nboot,replace=TRUE),nrow=nboot)
    bootdatay <- matrix(sample(ny,size=ny*nboot,replace=TRUE),nrow=nboot)
    bvec = vector(mode = "numeric", length = nboot)
    for(Bi in 1:nboot){
      mb = outer(x[bootdatax[Bi,]],y[bootdatay[Bi,]],"-")
      bvec[Bi] = hd(mb,q) + hd(mb,q2)
    }
    # p value
    p = mean(bvec>0)+.5*mean(bvec==0)
    p = 2*min(c(p,1-p))
    output[i,8] = p
    # confidence interval
    sbv = sort(bvec)
    ci = sbv[ilow]
    ci[2] = sbv[ihi]
    output[i,5] = ci[1]
    output[i,6] = ci[2]
  }
  temp = order(output[,8],decreasing=TRUE)
  zvec = alpha/c(1:length(qseq))
  output[temp,7] = zvec
  tmp <- data.frame(output)
  # add sig column
  # output$signif = rep("YES",nrow(output))
  # for(i in 1:nrow(output)){
  #   if(output[temp[i],8] > output[temp[i],7])output$signif[temp[i]]="NO"
  #   if(output[temp[i],8] <= output[temp[i],7])break
  # }
  out[[comp]] <- tmp
  names(out)[comp] <- paste0(gr_name1, " - ",gr_name2)
  }
  out
}

# ==========================================================================
#' Difference asymmetry function for two dependent groups or a single distribution
#'
#' \code{asymdhd} computes a difference asymmetry function for one or more
#' distributions. The difference asymmetry function provides perspective on the
#' degree a distribution is symmetric about zero, by quantifying the sum of q
#' and 1-q quantiles. If the groups do not differ, then
#' the difference scores should be symmetric about zero. In particular, the sum
#' of q and 1-q quantiles should be zero. If the distribution is symmetric the
#' function should be approximately a horizontal line. If in addition the median
#' of the difference scores is zero, the horizontal line will intersect the
#' y-axis at zero. Confidence intervals and p values are returned for each
#' quantile sum. The FWE is controlled via Hochberg's method, which is used to
#' determine critical p values based on the argument alpha. To plot the results
#' use \code{\link{plot_diff_asym}}.
#'
#' @param data A data frame in long format. One column is a factor describing the groups;
#'   another column contains the values/observations for each group. A properly formatted data
#'   frame can be created using \code{\link{mkt1}} or \code{\link{mkt2}}. Missing values are not
#'   allowed.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param q Quantiles to estimate - default = 0.05:0.05:.4 - must be <.5.
#' @param nboot Number of bootstrap samples - default = 1000.
#' @param alpha Expected long-run type I error rate - default = 0.05.
#' @param todo A list of comparisons to perform - default = NULL.
#' @param doall Set to TRUE to compute all comparisons - default = FALSE. Not
#'   executed if a \code{todo} list is provided.
#' @return A list of data frames, one data frame per comparison. Each data frame
#'   has one row per decile. The columns are: \itemize{
#'   \item Column 1 = quantiles
#'   \item Column 2 = quantiles of differences
#'   \item Column 3 = 1 - quantiles of differences
#'   \item Column 4 = sum of quantiles
#'   \item Column 5 = lower bounds of the confidence intervals
#'   \item Column 6 = upper bounds of the confidence intervals
#'   \item Column 7 = critical p_values based on Hochberg's method
#'   \item Column 8 = p_values
#'   }
#'
#' @section Note:
#' This function combines Rand Wilcox's difQpci and Dqdif R functions,
#' from Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/},
#' http://dornsife.usc.edu/labs/rwilcox/software/.
#'
#' @references
#' Wilcox, R.R. & Erceg-Hurn, D.M. (2012)
#' Comparing two dependent groups via quantiles.
#' J Appl Stat, 39, 2655-2664.
#'
#' @seealso \code{\link{hd}}
#' \code{\link{asymhd}} for independent groups
#'
#' @examples
#' set.seed(21) # generate data
#' n <- 100 # sample size per condition
#' C1 <- rnorm(n) # condition 1
#' C2 <- C1 + rnorm(n) + 2 # condition 2
#' # Data with 3 independent groups and 2 dependent conditions per group
#' library(tibble)
#' df <- tibble(gr = factor(c(rep("group1",n),rep("group2",n),rep("group3",n))),
#'   cond1 = c(C1, C1+rnorm(n), C1+rnorm(n)),
#'   cond2 = c(C2, C2 + 1, C2 + 3) ) # make tibble
#' library(dplyr)
#' df <- mutate(df, diff = cond1 - cond2)
#'
#' out <- asymdhd(df, diff ~ gr) # use the default parameters
#' out <- asymdhd(df, diff ~ gr, alpha = .90) # specify alpha
#' out <- asymdhd(df, diff ~ gr, nboot = 500) # specify the number of bootstrap samples
#' out <- asymdhd(df, diff ~ gr, todo = list("group1", "group2")) # specify list of comparisons
#' out <- asymdhd(df, diff ~ gr, q = seq(.1,.4,.1)) # specify the quantiles
#' out <- asymdhd(df, diff ~ gr, doall = TRUE) # compute all tests
#'
#' @export
asymdhd <- function(data = df,
                    formula = obs ~ gr,
                    qseq = seq(5,40,5)/100,
                    alpha = .05,
                    nboot = 1000,
                    todo = NULL,
                    doall = FALSE){
  # check input is a data frame
  if(!is.data.frame(data)){
    stop("data must be a data frame")
  }
  # subset data
  subf <- subset_formula(data, formula)
  if (length(todo)==0) { # no comparison is specified
    if (doall == FALSE) { # do not perform all tests
      if (length(subf$gr_names) > 1) {
        warning(paste0("Parameter column ",subf$param_col_name," contains more than 1 level. The difference asymmetry function is computed using observations from the first level only: ",subf$gr_names[1]))
      }
      todo <- list(subf$gr_names[1])
    }
    if (doall == TRUE) { # use all groups
      todo <- subf$gr_names
    }
  }
  # confidence interval's boundaries
  low <- round((alpha/2)*nboot)+1
  up <- nboot-low
  out <- vector("list", length(todo)) # declare list of difference asymmetry functions
  for(comp in 1:length(todo)){ # for each comparison
    dif <- data[data[[subf$param_col_name]] == todo[[comp]][1], subf$obs_col_name][[1]]
    nv = length(dif)
    gr_name1 <- todo[[comp]][1]
    output = matrix(0, nrow=length(qseq), ncol=8)
    dimnames(output) = list(NULL,c("quantile","Est_q","Est_1.minus.q","SUM","ci.low","ci.up","p_crit","p-value"))
    for(i in 1:length(qseq)){
      q <- qseq[i]
      output[i,1] <- q
      output[i,2] <- hd(dif,q=q)
      output[i,3] <- hd(dif,q=1-q)
      # output[i,4] = output[i,2] + output[i,3]
      # bootstrap samples
      bootdata <- matrix(sample(nv,size=nv*nboot,replace=TRUE),nrow=nboot)
      bvec <- vector(mode = "numeric", length = nboot)
      for(ib in 1:nboot){
        bvec[ib] <- hd(dif[bootdata[ib,]],q=q) + hd(dif[bootdata[ib,]],q=1-q)
      }
      # p value
      pv = mean(bvec<0) + .5*mean(bvec==0)
      p = 2*min(c(pv,1-pv))
      output[i,8] = p
      # confidence interval
      sbvec = sort(bvec)
      ci = sbvec[low]
      ci[2] = sbvec[up]
      output[i,5] = ci[1]
      output[i,6] = ci[2]
    }
    temp = order(output[,8],decreasing=TRUE)
    zvec = alpha/c(1:length(qseq))
    output[temp,7] = zvec
    output[,4] = output[,2] + output[,3]
    tmp <- data.frame(output)
    # add sig column
    # output$signif = rep("YES",nrow(output))
    # for(i in 1:nrow(output)){
    #   if(output[temp[i],8] > output[temp[i],7]) output$signif[temp[i]]="NO"
    #   if(output[temp[i],8] <= output[temp[i],7]) break
    # }
    out[[comp]] <- tmp
    names(out)[comp] <- gr_name1
  }
  out
}

