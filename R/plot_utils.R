#' Add bars marking quantiles to ggplot object created by \code{\link{plot_scat2}}.
#' Used in the README.md file to demonstrate the shift function.
#'
#' @param p A ggplot object returned by \code{\link{plot_scat2}}.
#' @param q_seq A sequence of quantiles - default = deciles.
#' @param col Colour of the bars
#' @param width Length of the bars
#' @param q_size Thickness of the bars, except for the median
#' @param md_size Thickness of the median bar
#' @param alpha Alpha transparency
#'
#' @return A ggplot object
#'
#' @export
plot_hd_bars <- function(p,
                         q_seq = seq(.1,.9,.1),
                         col = "grey21",
                         width = 0.5,
                         q_size = 0.5,
                         md_size = 1,
                         alpha = 1){
   nq <- length(q_seq)
   size_seq <- c(rep(q_size,floor(nq/2)), md_size, rep(q_size,floor(nq/2)))
   for (qi in 1:length(q_seq)){
    p <- p + geom_errorbar(stat = "summary",
                           fun.y = hd,
                           fun.ymin = hd,
                           fun.ymax = hd,
                           fun.args = list(q = q_seq[[qi]]),
                           colour = col,
                           width = width,
                           size = size_seq[[qi]],
                           alpha = alpha)
  }
  p
}

#' Add bar marking the mean
#' @export
plot_mean_bar <- function(p,
                          col = "grey21",
                          width = 0.5,
                          size = 0.5){
    p <- p + geom_errorbar(stat = "summary",
                           fun.y = mean,
                           fun.ymin = mean,
                           fun.ymax = mean,
                           colour = col,
                           width = width,
                           size = size)
  p
}



# =================================================================================
#' Add bars marking quantiles + links between quantiles to ggplot object created
#' by \code{\link{plot_scat2}}. Used in README to demonstrate the shift
#' function.
#'
#' @param p A ggplot object returned by \code{\link{plot_scat2}}.
#' Used in README to demonstrate the shift function.
#' @param sf A shift function for the 2 groups illustrated in p.
#' @param q_col Colour of the bars
#' @param q_width Length of the bars
#' @param q_size Thickness of the bars, except for the median
#' @param md_size Thickness of the median bar
#' @param link_col Colour of the links between quantiles - default "darkviolet"
#'   for negative differences, "darkorange2" for postivie differences
#' @param link_alpha Alpha transparency of the links between quantiles - default
#'   c(0.4, 1), maximum value for the median, decreasing towards the lowest
#'   value for the extreme quantiles.
#' @param add_rect Add rectangle marking the location of the quantiles in group 1.
#' @param rect_alpha Alpha transparency for the rectangle.
#' @param rect_col Colour for the rectangle
#' @param add_lab If TRUE, add labels for differences between extreme quantiles - default = FALSE
#' @param labres Number of decimales for the labels - default = 2
#' @param text_size Size of the labels - default = 5
#'
#' @export
plot_hd_links <- function(p,
                          sf = sf,
                          q_col = "grey21",
                          q_width = 0.5,
                          q_size = 0.5,
                          md_size = 1,
                          link_col = c("darkviolet","darkorange2"),
                          link_alpha = c(0.4, 1),
                          add_rect = FALSE,
                          rect_alpha = NULL,
                          rect_col = NULL,
                          add_lab = FALSE,
                          labres = 2,
                          text_size = 5){
  p <- plot_hd_bars(p,
                    col = q_col,
                    width = q_width,
                    q_size = q_size,
                    md_size = md_size,
                    q_seq = sf[[1]],
                    alpha = 1)
  # extract vectors from shift function -------------------------
  g1 <- sf[[2]] # group 1 deciles
  g2 <- sf[[3]] # group 2 deciles
  diff <- sf[[4]] # differences
  # create new variables ----------------------------------------
  diff_sign <- (sign(diff) > 0) + 1 # difference signs c(-1,1) -> c(1,2)
  q_seq <- sf[[1]]
  qn <- length(q_seq)
  deco <- c(seq(1,floor(qn/2)+1),seq(floor(qn/2),1)) # code of deciles
  alpha_seq <- seq(link_alpha[1], link_alpha[2], length.out = floor(qn/2)+1)
  line_size <- c(rep(q_size, floor(qn/2)),
                 md_size,
                 rep(q_size, floor(qn/2)))
  # add links ---------------------------------------------------
  for (d in 1:qn){
    p <- p + annotate("segment",
      x = 1 + q_width / 2,
      xend = 2 - q_width / 2,
      y = g1[d],
      yend = g2[d],
      colour = link_col[diff_sign[d]],
      alpha = alpha_seq[deco[d]],
      size = line_size[d])
  }
  # add rectangle
  if (add_rect == TRUE) {
    if (is.null(rect_alpha)){
      rect_alpha <- 0.2
    }
    if (is.null(rect_col)){
      rect_col <- "grey30"
    }
    p <- p + annotate("rect", xmin = 0.4, xmax = 1.25, ymin = g1[1], ymax = g1[qn],
      alpha = rect_alpha) # add rectangle
  }
  # add labels for differences between extreme deciles
  if (add_lab == TRUE){
    for (d in seq(1,qn,qn-1)){
      # if (diff_sign[d] == 1){
      #   cc <- link_col[2]
      # } else {
      #   cc <- link_col[1]
      # }
      p <- p + annotate("label",
        x = 1.5,
        y = min(g1[d],g2[d]) + abs(g1[d] - g2[d]) / 2,
        label = round(diff[d],labres),
        fill = link_col[diff_sign[d]],
        colour = "white",
        fontface = "bold",
        alpha = alpha_seq[deco[d]])
    } # for loop
  } # if add_lab
  p
}
# =================================================================================

#' Add difference labels to shift function
#'
#' Add labels of quantile difference values to one or more shift function plots.
#' Used in the README.md file to illustrate the shift function. Assumes an odd
#' number of quantiles with the median in the middle.
#'
#' @param p A list of ggplot objects generated by \code{\link{plot_sf}} or \code{\link{plot_pbsf}}.
#' @param sf A list of data frames generated by \code{\link{shifthd}},
#'   \code{\link{shiftdhd}}, \code{\link{shifthd_pbci}} or
#'   \code{\link{shiftdhd_pbci}}.
#' @param labres Number of decimales for the labels - default = 2.
#' @param link_col Label colours for negative and positive values.
#' @param link_alpha Alpha transparency of the labels - default = continuum between 0.4 and 1.
#' @param y_lab_nudge Amount by which to nudge the labels along the y axis.
#' @param text_size Text size - default = 5
#'
#' @export
add_sf_lab <- function(p,
                       sf,
                       labres = 2,
                       link_col = c("darkviolet","darkorange2"),
                       link_alpha = c(0.4, 1),
                       y_lab_nudge = 0,
                       text_size = 5){
  # check p and sf have same length
  if(length(p)!=length(sf)){
    stop("p and sf must have the same length")
  }
  # check p input is a list
  if(!is.list(p)){
    stop("p must be a list")
  }
  # check sf input is a list of data frames
  if(!is.list(sf)){
    stop("sf must be a list")
  }
  for (pc in 1:length(sf)) {
    if(!is.data.frame(sf[[pc]])){
      stop("input sf list must contain data.frames")
    }
  }
  plist <- vector("list", length(p)) # declare list of plot objects
  for (pc in 1:length(sf)) {
    # extract vectors from shift function -------------------------
    if(names(sf[[pc]][1]) == "q"){ # pbci shift function
      dec <- sf[[pc]][[1]]
      g1 <- sf[[pc]][[2]] # group 1 quantiles
      g2 <- sf[[pc]][[3]] # group 2 quantiles
      nq <- length(g1)
      deco <- c(seq(1,floor(nq/2)+1),seq(floor(nq/2),1)) # code of deciles
      alpha_seq <- seq(link_alpha[1], link_alpha[2], length.out = floor(nq/2)+1)
    } else { # pbse shift function
      dec <- seq(1,9,1)/10
      g1 <- sf[[pc]][[1]] # group 1 deciles
      g2 <- sf[[pc]][[2]] # group 2 deciles
      deco <- c(seq(1,5),seq(4,1)) # code of deciles
      alpha_seq <- seq(link_alpha[1], link_alpha[2], length.out = 5)
    }
    diff <- sf[[pc]]$difference # differences
    lo <- sf[[pc]]$ci_lower # lower confidence intervals
    hi <- sf[[pc]]$ci_upper # upper confidence intervals
    # create new variables ----------------------------------------
    diff_sign <- (sign(diff) > 0) + 1 # difference signs c(-1,1) -> c(1,2)
    # add labels for differences between extreme deciles
    for (d in 1:length(dec)){
      if (diff[d] < 0){ # negative difference
        nudge_y <- hi[d] + y_lab_nudge
      } else { # positive difference
        nudge_y <- lo[d] - y_lab_nudge
      }
      p[[pc]] <- p[[pc]] + annotate("label",
        x = g1[d],
        y = nudge_y,
        label = round(diff[d],labres),
        fill = link_col[diff_sign[d]],
        colour = "white",
        fontface = "bold",
        alpha = alpha_seq[deco[d]],
        size = text_size)
    } # quantile loop
  } # loop for list of ggplot objects
  p
}



