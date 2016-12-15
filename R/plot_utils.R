#' Quantile function to use with ggplot2 stat_summary
#'
plot_hd <- function(x,q){
  m <- hd(x,q)
  c(y = m, ymin = m, ymax = m)
}

plot_mean <- function(x){
  m <- mean(x)
  c(y = m, ymin = m, ymax = m)
}

#' Add bars marking quantiles to plot
#' Default to deciles.
#' To plot the deciles with a different bar size for the median,
#' use \code{\link{plot_dec_bars}}.
plot_hd_bars <- function(p,
                            q_seq = seq(.1,.9,.1),
                            col = "grey21",
                            width = 0.5,
                            size = 0.5){
  for (qi in 1:length(q_seq)){
    p <- p + geom_errorbar(stat = "summary",
                           fun.data = "plot_hd",
                           fun.args = list(q = q_seq(qi)),
                           colour = col,
                           width = width,
                           size = size)
  }
  p
}

#' Add bars marking the deciles to plot
#'
plot_dec_bars <- function(p,
                            q_seq = seq(.1,.9,.1),
                            col = "grey21",
                            width = 0.5,
                            dec_size = 0.5,
                            md_size = 1){
  size_seq <- c(rep(dec_size,4), md_size, rep(dec_size,4))
  for (qi in 1:length(q_seq)){
    p <- p + geom_errorbar(stat = "summary",
                           fun.data = "plot_hd",
                           fun.args = list(q = q_seq[[qi]]),
                           colour = col,
                           width = width,
                           size = size_seq[[qi]])
  }
  p
}

#' Add bars marking the quartiles to plot
#'
plot_quartile_bars <- function(p,
                               q_seq = c(.25, .5, .75),
                               col = "grey21",
                               width = 0.5,
                               q_size = 0.5,
                               md_size = 1){
  size_seq <- c(q_size, md_size, q_size)
  for (qi in 1:length(q_seq)){
    p <- p + geom_errorbar(stat = "summary",
                           fun.data = "plot_hd",
                           fun.args = list(q = q_seq[[qi]]),
                           colour = col,
                           width = width,
                           size = size_seq[[qi]])
  }
  p
}

#' Add bar marking the mean
#'
plot_mean_bar <- function(p,
                          col = "grey21",
                          width = 0.5,
                          size = 0.5){
    p <- p + geom_errorbar(stat = "summary",
                           fun.data = "plot_mean",
                           colour = col,
                           width = width,
                           size = size)
  p
}

# =================================================================================
#' Add deciles to plot + links between deciles
#'
#' Used in README to demonstrate the shift function.
#'
#' @export
plot_dec_links <- function(p,
                           sf = sf,
                           dec_col = "grey21",
                           dec_width = 0.5,
                           dec_size = 0.5,
                           md_size = 1,
                           link_col = c("darkviolet","darkorange2"),
                           link_alpha = c(0.4, 1),
                           add_rect = FALSE,
                           rect_alpha = NULL,
                           rect_col = NULL,
                           add_lab = FALSE,
                           labres = 2){
  p <- plot_dec_bars(p,
                     col = dec_col,
                     width = dec_width,
                     dec_size = dec_size,
                     md_size = md_size)
  # extract vectors from shift function -------------------------
  g1 <- sf[[1]] # group 1 deciles
  g2 <- sf[[2]] # group 2 deciles
  diff <- sf[[3]] # differences
  # create new variables ----------------------------------------
  diff_sign <- (sign(diff) > 0) + 1 # difference signs c(-1,1) -> c(1,2)
  deco <- c(seq(1,5),seq(4,1)) # code of deciles
  alpha_seq <- seq(link_alpha[1], link_alpha[2], length.out = 5)
  line_size <- c(dec_size,dec_size,dec_size,dec_size,
                 md_size,
                 dec_size,dec_size,dec_size,dec_size)
  # add links ---------------------------------------------------
  for (d in 1:9){
  p <- p + annotate("segment",
                    x = 1 + dec_width / 2,
                    xend = 2 - dec_width / 2,
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
  p <- p + annotate("rect", xmin = 0.4, xmax = 1.25, ymin = g1[1], ymax = g1[9],
                       alpha = rect_alpha) # add rectangle
  }
  # add labels for differences between extreme deciles
  if (add_lab == TRUE){
    for (d in seq(1,9,8)){
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
#' Add labels of difference values to shift function plot.
#' Used in README.md file to illustrate the shift function.
#'
#' @export
add_sf_lab <- function(p,
                       sf,
                       dec = seq(1,9),
                       labres = 2,
                       link_col = c("darkviolet","darkorange2"),
                       link_alpha = c(0.4, 1),
                       y_lab_nudge = 0){
  # extract vectors from shift function -------------------------
  g1 <- sf[[1]] # group 1 deciles
  g2 <- sf[[2]] # group 2 deciles
  diff <- sf[[3]] # differences
  lo <- sf[[4]] # lower confidence intervals
  hi <- sf[[5]] # upper confidence intervals
  # create new variables ----------------------------------------
  diff_sign <- (sign(diff) > 0) + 1 # difference signs c(-1,1) -> c(1,2)
  deco <- c(seq(1,5),seq(4,1)) # code of deciles
  alpha_seq <- seq(link_alpha[1], link_alpha[2], length.out = 5)
  # add labels for differences between extreme deciles
    for (d in 1:length(dec)){
      if (diff[dec[d]] < 0){ # negative difference
        nudge_y <- hi[dec[d]] + y_lab_nudge
      } else { # positive difference
        nudge_y <- lo[dec[d]] - y_lab_nudge
      }
      p <- p + annotate("label",
                        x = g1[dec[d]],
                        y = nudge_y,
                        label = round(diff[dec[d]],labres),
                        fill = link_col[diff_sign[dec[d]]],
                        colour = "white",
                        fontface = "bold",
                        alpha = alpha_seq[deco[dec[d]]])
    } # for loop
  p
}

#' Add text annotation of quartiles
#'
annotate_quartiles <- function(p,
                               data,
                               x = 0,
                               hjust = 0,
                               vjust = 0,
                               size = 10){
  hdq = c(hd(out,.25),hd(out,.5),hd(out,.75)) # compute quartiles
  caption <- as.character(round(hdq, digits=2)) # turn into characters
  p <- p + annotate("text", x = x, y = hdq[1], label = caption[1],
                    hjust = hjust, vjust = vjust, size = size) +
    annotate("text", x = x, y = hdq[2], label = caption[2],
             hjust = hjust, vjust = vjust, size = size) +
    annotate("text", x = x, y = hdq[3], label = caption[3],
             hjust = hjust, vjust = vjust, size = size)
  p
}

