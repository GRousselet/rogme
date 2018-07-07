#' Plot one-dimensional scatterplots for 2 groups
#'
#' \code{plot_scat2} produces scatterplots for 2 marginal distributions.
#' The scatterplots are jittered using \code{\link[ggbeeswarm]{geom_quasirandom}}.
#'
#' @param data A data frame in long format. One column is a factor describing the groups;
#'   another column contains the values/observations for each group. A properly formatted data
#'   frame can be created using \code{\link{mkt2}}. Missing values are not
#'   allowed.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param xlabel Option to set different name - default NULL to use data frame column names.
#' @param ylabel Option to set different name - default NULL to use data frame column names.
#' @param ... Input arguments for ggbeeswarm::geom_quasirandom
#'
#' @return A ggplot object.
#'
#' @examples
#' # generate data
#' set.seed(21)
#' g1 <- rnorm(1000) + 6
#' g2 <- rnorm(1000) * 1.5 + 6
#'
#' # make tibble
#' df <- mkt2(g1, g2)
#' # make scatterplots
#' ps <- plot_scat2(data = df,
#'   formula = obs ~ gr,
#'   xlabel = "",
#'   ylabel = "Scores (a.u.)",
#'   alpha = 1,
#'   shape = 21,
#'   colour = "grey10",
#'   fill = "grey90") # scatterplots
#' ps <- ps + coord_flip()
#' ps
#'
#' @export
plot_scat2 <- function(data = df,
                       formula = obs ~ gr,
                       xlabel = NULL,
                       ylabel = NULL,
                       ...){
  # subset formula
  subf <- subset_formula(data, formula)
  xplot = subf$param_col_name
  yplot = subf$obs_col_name
  p <- ggplot(data, aes_string(x = xplot, y = yplot, fill = xplot,
                               colour = xplot, shape = xplot))
  p <- p + ggbeeswarm::geom_quasirandom(...) +
    theme_bw() +
    # scale_colour_manual(values = symb_col) +
    # scale_fill_manual(values = symb_fil) +
    # scale_shape_manual(values = symb_shape) +
    theme(legend.position = "none") +
    theme(axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold"),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14))
  # override axis labels
  if (!is.null(xlabel)){
    p <- p + xlab(xlabel)
  }
  if (!is.null(ylabel)){
    p <- p + ylab(ylabel)
  }
  p
}

#' Plot quantiles and confidence intervals
#'
#' Using the output of \code{\link{quantiles_pbci}}, create a ggplot object
#' showing specified quantiles (default to deciles) and their 95% percentile
#' bootstrap confidence intervals. A vertical line marks the median.
#'
#' @param data Data frame created by `quantiles_pbci`.
#' @param qseq Sequence of quantiles to plot - assumes median is in the middle of the sequence - default = deciles.
#'
#' @seealso \code{\link{quantiles_pbci}}
#'
#' @examples
#' set.seed(7)
#' # make fake skewed data
#' x <- rgamma(100, shape=3, scale=1)*10+250
#' # compute quantiles and their percentile bootstrap confidence intervals
#' out <- quantiles_pbci(x,q=seq(1,9)/10,nboot=2000,alpha=0.05)
#' # make decile plot
#' p <- plot_hd_ci(data=out,plotzero=TRUE,label.x="Onsets in ms",
#'                    hjust=-.05,vjust=.5,size_text=5,
#'                    colour_dec = "grey10",fill_dec = "grey90",
#'                    colour_line = "grey10", linetype_line = 1, size_line = 1) +
#'                    scale_y_continuous(limits=c(250, 350),breaks=seq(250,350,25))
# p
#'
#' @export
plot_hd_ci <- function(data = out,
                       qseq = seq(.1,.9,.1),
                        plotzero = TRUE,
                        label.x = "Differences",
                        hjust = -.05,
                        vjust = .2,
                        size_text = 6,
                        colour_q = "#009E73",
                        fill_q = "white",
                        size_q = 4,
                        shape_q = 21,
                        colour_line = "#009E73",
                        size_line = 1,
                        linetype_line = 2,
                        colour_zero = "black",
                        size_zero = .5,
                        linetype_zero = 1){
  md.loc <- floor(length(data$quantile)/2) + 1
  md <- data$est_q[md.loc] # median
  md.c <- as.character(round(md, digits=1)) # turn into characters
  lo.c <- as.character(round(data$ci.low[md.loc], digits=1)) # turn into characters
  up.c <- as.character(round(data$ci.up[md.loc], digits=1)) # turn into characters
  caption <- paste("Median = \n ",md.c," [",lo.c,", ",up.c,"]",sep="")
  if (all.equal(qseq,seq(.1,.9,.1))){
    label.y <- "Deciles"
  }  else {
    label.y <- "Quantiles"
  }
  p <- ggplot(data=out, aes(x=quantile*10, y=est_q)) +
    geom_abline(intercept = md, slope = 0,
                colour = colour_line,
                size = size_line,
                linetype = linetype_line)
  if (plotzero){
    p <- p + geom_abline(intercept = 0, slope = 0,
                         colour = colour_zero,
                         size = size_zero,
                         linetype = linetype_zero)
  }
  p <- p + geom_linerange(aes(ymin=ci.low, ymax=ci.up), colour=colour_line,size=1) +
    geom_point(colour = colour_q,
               size = size_q,
               shape = shape_q,
               fill = fill_q) +
    theme_bw() +
    labs(x=label.y) +
    labs(y=label.x) +
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold")) +
    scale_x_continuous(breaks=seq(1,9,1)) +
    annotate("text", x = 5, y = data$ci.up[5], label = caption[1],
             hjust = hjust, vjust = vjust, size = size_text) +
    coord_flip()
  p
}





# ----------------------------------------------------------------------------
#' Plot paired observations
#'
#' Scatterplot of paired observations with reference line of no effect.
#' Quartiles of each condition are superimposed. Quartiles are estimated using the
#' Harrell-Davis estimator.
#' @param df Data frame with paired observations in two columns.
#' @param formula A formula with format response variable ∼ predictor variable,
#'   where ~ (tilde) means "is modeled as a function of".
#' @param axis.steps Steps between x and y tick marks - default = 2.
#' @param min.x,min.y,max.x,max.y Specify axis limits - default square axes
#' @param colour_p Colour parameter of the scatterplot - default "black".
#' @param size_p Size parameter of the scatterplot - default = 5.
#' @param stroke_p Stroke parameter of the scatterplot - default = 1,
#' @param shape_p Shape parameter of the scatterplot - default = 21,
#' @param colour_p Colour parameter of the scatterplot - default = "black",
#' @param fill_p Fill parameter of the scatterplot - default = "#ffb347",
#' @param alpha_p Alpha parameter of the scatterplot - default = .5,
#' @param linetype_q Linetype of the segments marking the quartiles - default = "dashed",
#' @param size_q Size of the segments marking the quartiles - default = 1,
#' @param alpha_q Alpha of the segments marking the quartiles - default = .5,
#' @param colour_q Colour of the segments marking the quartiles - default = "black"
#' @examples
#' df <- tibble(cond1 = rnorm(50), cond2 = cond1 + rnorm(50))
#' plot_scat2d(df, formula = cond2 ~ cond1) # basic call
#' plot_scat2d(df, formula = cond2 ~ cond1, size_q=3) # specify size of quartile segments
#' plot_scat2d(df, formula = cond2 ~ cond1, size_q=c(1,2,1)) # use thicker line for median
#' plot_scat2d(df, formula = cond2 ~ cond1, linetype_q = "longdash") # specify linetype - default = dashed
#' @seealso \code{\link{hd}}
#' @export
plot_scat2d <- function(df = df,
                        formula = cond2 ~ cond1,
                        min.x = NA,
                        min.y = NA,
                        max.x = NA,
                        max.y = NA,
                        axis.steps = 2,
                        size_p = 5,
                        stroke_p = 1,
                        shape_p = 21,
                        colour_p = "black",
                        fill_p = "#ffb347",
                        alpha_p = .5,
                        linetype_q = "dashed",
                        size_q = 1,
                        alpha_q = .5,
                        colour_q = "black"){
  # subset formula
  subf <- subset_formula_wide(df, formula)
  xplot = subf$x_col_name
  yplot = subf$y_col_name
  # make data.frames for plotting quartile segments
  hd1.25 <- hd(df[[xplot]],.25)
  hd1.5 <- hd(df[[xplot]],.5)
  hd1.75 <- hd(df[[xplot]],.75)
  hd2.25 <- hd(df[[yplot]],.25)
  hd2.5 <- hd(df[[yplot]],.5)
  hd2.75 <- hd(df[[yplot]],.75)
  df.25 <- data.frame(hd1=hd1.25,hd2=hd2.25)
  df.5 <- data.frame(hd1=hd1.5,hd2=hd2.5)
  df.75 <- data.frame(hd1=hd1.75,hd2=hd2.75)

  # quartile plot parameters
  if(length(linetype_q)==1){
    linetype_q = rep(linetype_q,3)
  }
  if(length(size_q)==1){
    size_q = rep(size_q,3)
  }
  if(length(alpha_q)==1){
    alpha_q = rep(alpha_q,3)
  }
  if(length(colour_q)==1){
    colour_q = rep(colour_q,3)
  }

  # plot limits
  if (is.na(min.x)){
    min.x <- min(df[[xplot]],df[[yplot]])
  }
  if (is.na(max.x)){
    max.x <- max(df[[xplot]],df[[yplot]])
  }
  if (is.na(min.y)){
    min.y <- min(df[[xplot]],df[[yplot]])
  }
  if (is.na(max.y)){
    max.y <- max(df[[xplot]],df[[yplot]])
  }

  # scatterplot of paired observations -----------------
  p <- ggplot(df, aes_string(x = xplot, y = yplot)) +
    # reference line
    geom_abline(intercept = 0) +
    # quartiles
    scale_x_continuous(breaks=seq(floor(min.x),ceiling(max.x),axis.steps)) +
    scale_y_continuous(breaks=seq(floor(min.y),ceiling(max.y),axis.steps)) +
    coord_cartesian(xlim = c(floor(min.x), ceiling(max.x)),
                    ylim = c(floor(min.y), ceiling(max.y))) +
    geom_segment(aes(x=hd1,y=min.y*2,xend=hd1,yend=hd2),data=df.25,linetype=linetype_q[1],size=size_q[1],alpha=alpha_q[1],colour=colour_q[1]) +
    geom_segment(aes(x=hd1,y=min.y*2,xend=hd1,yend=hd2),data=df.5,linetype=linetype_q[2],size=size_q[2],alpha=alpha_q[2],colour=colour_q[2]) +
    geom_segment(aes(x=hd1,y=min.y*2,xend=hd1,yend=hd2),data=df.75,linetype=linetype_q[3],size=size_q[3],alpha=alpha_q[3],colour=colour_q[3]) +
    geom_segment(aes(x=min.x*2,y=hd2,xend=hd1,yend=hd2),data=df.25,linetype=linetype_q[1],size=size_q[1],alpha=alpha_q[1],colour=colour_q[1]) +
    geom_segment(aes(x=min.x*2,y=hd2,xend=hd1,yend=hd2),data=df.5,linetype=linetype_q[2],size=size_q[2],alpha=alpha_q[2],colour=colour_q[2]) +
    geom_segment(aes(x=min.x*2,y=hd2,xend=hd1,yend=hd2),data=df.75,linetype=linetype_q[3],size=size_q[3],alpha=alpha_q[3],colour=colour_q[3]) +
    # scatterplot
    geom_point(size=size_p,
               stroke=stroke_p,
               shape=shape_p,
               colour=colour_p,
               fill=fill_p,
               alpha=alpha_p) +
    theme_bw() +
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold"),
          legend.title = element_blank(),
          plot.title = element_text(size=20)) +
    labs(title="Paired observations")
  p
}
# ----------------------------------------------------------------------------
