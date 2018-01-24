#' Plot difference asymmetry function
#'
#' Plot difference asymmetry function using output from
#' Rand Wilcox's qwmwhd or difQpci functions. For naming consistency,
#' these functions have been renamed `asymhd` and `asymdhd` in `rogme`.
#' @export
plot_diff_asym <- function(data = df,
                           xlabel = "Quantiles",
                           ylabel = "Quantile sum = q + 1-q",
                           symb_col = NULL,
                           symb_fill = NULL,
                           symb_size = 5,
                           symb_shape = 21,
                           diffint_col = NULL,
                           diffint_size = .5,
                           dec_line_col = NULL,
                           dec_line_alpha = .5,
                           dec_line_size = 1.5){
  ylim <- max(max(abs(data$ci.up)),max(abs(data$ci.low)))
  ylim <- c(-ylim,ylim)
  data$deco <- seq(1,nrow(data))
  if (is.null(symb_col)){
    symb_col <- "black"
  }
  if (is.null(symb_fill)){
    symb_fill <- c("white","grey30")
  }
  if (is.null(diffint_col)){
    diffint_col <- "black"
  }
  if (is.null(dec_line_col)){
    dec_line_col <- "grey50"
  }

  p <- ggplot(data, aes_string(x="quantile", y="SUM")) +
    geom_hline(yintercept=0,linetype=2,alpha=0.5) + # x=0 reference line
    # vertical bars for uncertainty intervals
    geom_linerange(aes(ymin=ci.low, ymax=ci.up), colour = diffint_col,
                   size = diffint_size) +
    # line joining the quantiles
    geom_line(colour = dec_line_col, alpha = dec_line_alpha, linetype = "solid",
              size = dec_line_size) +
    # symbols marking the quantiles
    geom_point(aes(fill = deco), colour="black", size = symb_size, shape = symb_shape) +
    scale_fill_gradient(low = symb_fill[1], high = symb_fill[2], guide = FALSE) +
    xlab(xlabel) +
    ylab(ylabel) +
    theme_bw() +
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold")) +
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(breaks = data$quantile)
  p
}




