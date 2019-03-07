#' Hierarchical shift function for one group, two dependent conditions
#'
#' \code{hsf} returns a hierarchical shift function for one group of
#' participants, tested in two dependent conditions (see \href{https://github.com/GRousselet/rogme/blob/master/docs/hsf.md}{vignette} on github). Full distributions of
#' measurements must be available for each participant and condition. First,
#' quantiles are computed for the distribution of measurements from each condition and each
#' participant. Second, the quantiles are subtracted in each participant. Third,
#' a trimmed mean is computed across participants for each quantile. Confidence
#' intervals and p values are also computed. Correction for multiple comparisons
#' across quantiles is achieved using Hochberg's 1988 procedure. Plot the shift
#' function using \code{plot_hsf}.
#'
#' @references Rousselet, G. A., & Wilcox, R. R. (2019, January 17). Reaction
#'   times and other skewed distributions: problems with the mean and the
#'   median. https://doi.org/10.31234/osf.io/3y54r
#'
#'   Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of
#'   significance. Biometrika, 75(4), 800-802.
#'
#' @param data A data frame in long format. Missing values are not allowed.
#' @param formula A formula with format response variable ∼ predictor variable + id,
#'   where ~ (tilde) means "is modeled as a function of" and '+ id' indicates the variable containing the participants' id number.
#' @param qseq Quantiles to estimate - default = deciles.
#' @param tr Percentage of trimming, value between 0 and 1 - default = 0.2 = 20\%. Set to zero to get results for the mean.
#' @param alpha Alpha level - default 0.05.
#' @param qtype Type of quantile estimation algorithm to pass to \code{quantile} - default = 8.
#' @param todo Order of the groups to compare - default = 1 minus 2.
#' @param null.value Null value to compute P values for the quantile differences - default = 0.
#' @param adj_method Name of method to adjust for multiple quantile comparisons, passed to \code{p.adjust} - default = "hochberg".
#'
#' @return A list of 8 results:
#' \enumerate{
#'   \item \strong{comparison}: names of two conditions being compared.
#'   \item \strong{individual_sf}: shift functions for every participant.
#'   \item \strong{group_differences}: group quantile differences.
#'   \item \strong{ci}: group confidence intervals for quantile differences.
#'   \item \strong{pvalues}: P values for every difference.
#'   \item \strong{adjusted_pvalues}: P values adjusted for multiple comparisons.
#'   \item \strong{null_value}: null value used to compute P values.
#'   \item \strong{quantiles}: quantiles estimated in each participant and condition.
#'   }
#'
#' @seealso \code{\link{plot_hsf}} to plot the results.
#'
#' @examples
#' set.seed(22) # subset random sample of participants from the French Lexicon Project
#' id <- unique(flp$participant)
#' df <- subset(flp, flp$participant %in% sample(id, 50, replace = FALSE))
#' out <- hsf(df, rt ~ condition + participant) # use the default parameters
#' plot_hsf(out) # plot results. Shift functions are overall negative, as
#' participants tend to be faster in the Word condition than in the Non-Word
#' condition.
#'
#' out <- hsf(df, rt ~ condition + participant, qseq = c(.25, .5, .75)) # estimate quartiles only
#'
#' out <- hsf(df, rt ~ condition + participant, todo = c(2,1)) # reverse comparison
#'
#' @export
hsf <- function(data = df,
                formula = obs ~ cond + id,
                qseq = seq(0.1,0.9,0.1),
                tr = 0.2,
                alpha = 0.05,
                qtype = 8,
                todo = c(1,2),
                null.value = 0,
                adj_method = "hochberg"){
  # Check input is a data frame
  if(!is.data.frame(data)){
    stop("data must be a data frame")
  }
  # Subset data
  subf <- subset_formula_hsf(data, formula)
  if (length(subf$gr_names) > 2) {
    warning(paste0("Parameter column ",subf$param_col_name," contains more than 2 levels. The shift function is computed based on the first 2 levels: ",subf$gr_names[1], " vs. ",subf$gr_names[2]))
  }
  out <- vector("list", 8) # declare list of results

  # Change order of conditions?
  tmp <- subf$gr_names
  subf$gr_names[1] <- tmp[todo[1]]
  subf$gr_names[2] <- tmp[todo[2]]

  # Compute shift functions for all participants
  data.s <- subset(data, data[subf$param_col_name] == subf$gr_names[1])
  q.1 <- matrix(unlist(tapply(data.s[[subf$obs_col_name]], list(data.s[[subf$id_col_name]]), quantile, probs = qseq, type = qtype)), nrow=length(qseq))
  data.s <- subset(data, data[subf$param_col_name] == subf$gr_names[2])
  q.2 <- matrix(unlist(tapply(data.s[[subf$obs_col_name]], list(data.s[[subf$id_col_name]]), quantile, probs = qseq, type = qtype)), nrow=length(qseq))
  q.diff <- q.1 - q.2 # quantile differences

  # Compute group trimmed means
  group_diff = apply(q.diff, 1, mean, trim = tr)

  # Compute confidence intervals for group trimmed means
  group_ci <- apply(q.diff, 1, trimci, tr = tr, alpha = alpha) # 2 x length(qseq)

  # P values
  pval <- apply(q.diff, 1, trimpval, tr = tr, null.value = null.value) # length(qseq)

  # Adjusted P values
  adjust_pval <- p.adjust(pval, method = adj_method)

  names(out) <- c('comparison','individual_sf', 'group_differences', 'ci',
                  'pvalues', 'adjusted_pvalues', 'null_value', 'quantiles')
  out[[1]] <- paste0(subf$gr_names[1], " - ",subf$gr_names[2])
  out[[2]] <- q.diff
  out[[3]] <- group_diff
  out[[4]] <- group_ci
  out[[5]] <- pval
  out[[6]] <- adjust_pval
  out[[7]] <- null.value
  out[[8]] <- qseq
  out
}

#' Percentile bootstrap hierarchical shift function for one group, two dependent conditions
#'
#' \code{hsf_pb} returns a percentile bootstrap hierarchical shift function for
#' one group of participants, tested in two dependent conditions (see
#' \href{https://github.com/GRousselet/rogme/blob/master/docs/hsf.md}{vignette}
#' on github). Full distributions of measurements must be available for each
#' participant and condition. First, quantiles are computed for the distribution
#' of measurements from each condition and each participant. Second, the
#' quantiles are subtracted in each participant. Third, a trimmed mean is
#' computed across participants for each quantile. Confidence intervals are
#' computed using the percentile bootstrap. This function focuses on estimation, so
#' there is no correction for multiple comparisons and no p values. Plot the
#' shift function using \code{plot_hsf_pb} or \code{plot_hsf_pb_dist}.
#'
#' @references Rousselet, G. A., & Wilcox, R. R. (2019, January 17). Reaction
#'   times and other skewed distributions: problems with the mean and the
#'   median. https://doi.org/10.31234/osf.io/3y54r
#'
#' @param data A data frame in long format. Missing values are not allowed.
#' @param formula A formula with format response variable ∼ predictor variable + id,
#'   where ~ (tilde) means "is modeled as a function of" and '+ id' indicates the variable containing the participants' id number.
#' @param qseq Quantiles to estimate - default = deciles.
#' @param tr Percentage of trimming, value between 0 and 1 - default = 0.2 = 20\%. Set to zero to get results for the mean.
#' @param alpha Alpha level - default 0.05.
#' @param qtype Type of quantile estimation algorithm to pass to \code{quantile} - default = 8.
#' @param todo Order of the groups to compare - default = 1 minus 2.
#' @param nboot Number of boostrap samples - default = 1000.
#'
#' @return A list of 8 results:
#' \enumerate{
#'   \item \strong{comparison}: names of two conditions being compared.
#'   \item \strong{individual_sf}: shift functions for every participant.
#'   \item \strong{group_differences}: group quantile differences.
#'   \item \strong{ci}: group confidence intervals for quantile differences.
#'   \item \strong{hdi}: group highest density intervals for quantile differences.
#'   \item \strong{quantiles}: quantiles estimated in each participant and condition.
#'   \item \strong{boot_samples}: bootstrap differences for each quantile.
#'   \item \strong{nboot}: number of bootstrap samples.
#'   }
#'
#' @seealso \code{\link{plot_hsf_pb}} to plot the results.
#' \code{\link{plot_hsf_pb_dist}} to plot the distributions of bootstrap samples.
#'
#' @examples
#' set.seed(22) # subset random sample of participants from the French Lexicon Project
#' id <- unique(flp$participant)
#' df <- subset(flp, flp$participant %in% sample(id, 50, replace = FALSE))
#' out <- hsf_pb(df, rt ~ condition + participant) # use the default parameters
#' plot_hsf_pb(out) # plot results. Shift functions are overall negative, as
#' participants tend to be faster in the Word condition than in the Non-Word
#' condition.
#'
#' out <- hsf_pb(df, rt ~ condition + participant, qseq = c(.25, .5, .75)) # estimate quartiles only
#'
#' out <- hsf_pb(df, rt ~ condition + participant, todo = c(2,1)) # reverse comparison
#'
#' @export
hsf_pb <- function(data = df,
                    formula = obs ~ cond + id,
                    qseq = seq(0.1,0.9,0.1),
                    tr = 0.2,
                    alpha = 0.05,
                    qtype = 8,
                    todo = c(1,2),
                    nboot = 1000){
  # Check input is a data frame
  if(!is.data.frame(data)){
    stop("data must be a data frame")
  }
  # Subset data
  subf <- subset_formula_hsf(data, formula)
  if (length(subf$gr_names) > 2) {
    warning(paste0("Parameter column ",subf$param_col_name," contains more than 2 levels. The shift function is computed based on the first 2 levels: ",subf$gr_names[1], " vs. ",subf$gr_names[2]))
  }
  out <- vector("list", 8) # declare list of results
  
  # Define quantiles
  icrit <- round((1-alpha)*nboot) 
  ilo <- round((alpha/2)*nboot)
  iup <- nboot - ilo
  ilo <- ilo + 1
  
  # Change order of conditions?
  tmp <- subf$gr_names
  subf$gr_names[1] <- tmp[todo[1]]
  subf$gr_names[2] <- tmp[todo[2]]
  
  # simplify data frame - make tibble
  df <- tibble(id = data[[subf$id_col_name]],
               obs = data[[subf$obs_col_name]],
               cond = data[[subf$param_col_name]])
  
  # Number of participants
  np <- length(unique(df$id))
  
  #use filter to split conditions:
  #   data1 <- filter(data, conditions == cond1)
  # 
  # df <- tibble(rt = rnorm(33), id = factor(c(rep(1,10),rep(2,11),rep(3,12))))
  # df <- mutate(group_by(df, id), group_row = 1:n())
  # data.matrix(spread(df, id, rt))
  
  unique_id <- as.vector(unique(data[[subf$id_col_name]]), mode = "numeric")
 
  # Compute shift functions for all participants
  df1 <- dplyr::filter(df, df$cond == subf$gr_names[1])
  nt1 <- unlist(tapply(df1$obs, df1$id, length))
  q.1 <- matrix(unlist(tapply(df1$obs, df1$id, quantile, probs = qseq, type = qtype)), nrow=length(qseq))
  
  df2 <- dplyr::filter(df, df$cond == subf$gr_names[2])
  nt2 <- unlist(tapply(df2$obs, df2$id, length))
  q.2 <- matrix(unlist(tapply(df2$obs, df2$id, quantile, probs = qseq, type = qtype)), nrow=length(qseq))
  
  q.diff <- q.1 - q.2 # quantile differences
 
  # Compute group trimmed means
  group_diff = apply(q.diff, 1, mean, trim = tr)

  # Compute confidence intervals for group trimmed means ---------------------
  boot_qdiff <- array(data = 0, dim = c(nboot, length(qseq), np))
  
  # turn data frames into matrices for bootstrapping
  df1 <- dplyr::mutate(dplyr::group_by(df1, id), group_row = 1:dplyr::n())
  mat1 <- data.matrix(tidyr::spread(df1, id, obs))
  mat1 <- mat1[,-c(1,2)] # trials x participants
  
  df2 <- dplyr::mutate(dplyr::group_by(df2, id), group_row = 1:dplyr::n())
  mat2 <- data.matrix(tidyr::spread(df2, id, obs))
  mat2 <- mat2[,-c(1,2)] # trials x participants
  
  # bootstrap participants
  boot_id <- matrix(sample(unique_id, np * nboot, replace = TRUE), nrow = nboot)
  
  for(B in 1:nboot){
    
    for(CP in 1:np){ # bootstrap trials for each bootstrap participant
      boot_data1 <- sample(na.omit(mat1[,boot_id[B,CP]]), nt1[boot_id[B,CP]], replace = TRUE)
      
      boot_data2 <- sample(na.omit(mat2[,boot_id[B,CP]]), nt2[boot_id[B,CP]], replace = TRUE)
      
      boot_qdiff[B,,CP] <- quantile(boot_data1, probs = qseq, type = 8, names = FALSE) - quantile(boot_data2, probs = qseq, type = 8, names = FALSE)
    }
    
  }
  
  boot_tm <- apply(boot_qdiff, c(1,2), mean, trim = tr) # nboot x length(qseq)
  sort_boot_tm <- apply(boot_tm, 2, sort)
  boot_ci <- matrix(data = 0, nrow = 2, ncol = length(qseq))
  boot_ci[1,] <- sort_boot_tm[ilo,]
  boot_ci[2,] <- sort_boot_tm[iup,]
  
  boot_hdi <- apply(boot_tm, 2, HDInterval::hdi, credMass = 1-alpha)
  # -----------------------------------------------------------------

  names(out) <- c('comparison','individual_sf', 'group_differences', 'ci',
                  'hdi', 'quantiles', 'bootstrap_samples', 'nboot')
  out[[1]] <- paste0(subf$gr_names[1], " - ",subf$gr_names[2])
  out[[2]] <- q.diff
  out[[3]] <- group_diff
  out[[4]] <- boot_ci
  out[[5]] <- boot_hdi
  out[[6]] <- qseq
  out[[7]] <- boot_tm
  out[[8]] <- nboot
  out
}