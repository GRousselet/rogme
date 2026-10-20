#' @import ggplot2
#' @importFrom tibble tibble
#' @importFrom stats complete.cases na.omit p.adjust pbeta pt qnorm qt quantile var
#' @importFrom utils combn
NULL

# column names used in non-standard evaluation (tibble/ggplot2/tidybayes)
utils::globalVariables(c("boot_samp", "ci.low", "ci.up", "ci_lower", "ci_upper",
                         "dec", "deco", "est_q", "hd1", "hd2", "id", "obs",
                         "participants", "ymax", "ymin"))
