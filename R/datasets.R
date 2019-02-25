#' Reaction time data from the French Lexicon Project
#'
#' A dataset containing reaction times data from 959 participants in 2 conditions
#' from a lexical decision task. File created using the script getflprtdata.Rmd
#' available in /data-raw/.
#'
#' @format A data frame with 1917983 rows and 4 variables:
#' \describe{
#'   \item{participant}{id factor, from 1 to 959}
#'   \item{rt}{single trial reaction times, in ms}
#'   \item{acc}{accuracy, 1 for correct and 0 for incorrect responses}
#'   \item{condition}{condition factor, with levels "word" and "non-word"}
#' }
#' @source \url{https://sites.google.com/site/frenchlexicon/results}
"flp"

#' ERP onsets from 120 participants
#'
#' A dataset containing estimated onsets of face - noise ERP differences
#' from 120 participants
#'
#' @format A vector of 120 values in ms.
#'
#' @source \url{https://onlinelibrary.wiley.com/doi/full/10.1111/ejn.13100}
"onsets"

#'  Fake data from 35 participants
#'
#'  A dataset containing fake data from 35 participants in 2 conditions.
#'
#'  @format A data frame with 35 rows and 3 variables:
#'  \describe{
#'    \item{participant}{id factor, from 1 to 35}
#'    \item{condition1}{results from condition 1, in arbitrary units}
#'    \item{condition2}{results from condition 2, in arbitrary units}
#'  }
#'
#'  @source \url{https://www.biorxiv.org/content/10.1101/121079v2}
"pdata"
