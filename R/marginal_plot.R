#' Plot 2 groups: KDE + rug + deciles
#'
#' Plot kernel density estimates + rug plots + deciles
#' for 2 groups stored in a data frame.
#'
#' @export
plot_kde_rug_dec2 <- function(data = df){
  cdat <- plyr::ddply(data, "gr", summarise, deciles = q1469(obs))
  hd05 <- plyr::ddply(data, "gr", summarise, hd = hd(obs,0.5))
  #cc <- "grey80" # colour to plot deciles
  p <- ggplot(data, aes(x=obs, fill=gr)) + geom_density(alpha=.3) +
    facet_grid(gr ~ .) +
    geom_vline(data=hd05, aes(xintercept=hd,  colour=gr),
               linetype="solid", size=2, alpha=.5) + # thicker median
    geom_vline(data=cdat, aes(xintercept=deciles,  colour=gr),
               linetype="solid", size=1, alpha=.5) +
    geom_rug() +
    theme_bw() +
    theme(legend.position="none",
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold"),
          strip.text.y = element_text(size = 20, colour = "white"),
          strip.background = element_rect(colour="darkgrey", fill="darkgrey")) +
    ylab("Density")
  p
}

#' Plot 1 group: KDE + rug + deciles
#'
#' Plot kernel density estimate + rug plot + superimposed deciles
#' for 1 group stored in a data frame
#'
#' @export
plot_kde_rug_dec1 <- function(data=df,fill.colour="grey30",fill.alpha=.3){
  cdat <- plyr::ddply(data, "gr", summarise, deciles=q1469(obs))
  hd05 <- plyr::ddply(data, "gr", summarise, hd=hd(obs,0.5))
  cc <- "grey80" # colour to plot deciles
  p <- ggplot(data, aes(x=obs)) +
    geom_density(alpha=fill.alpha,fill=fill.colour,colour="black") +
    geom_vline(xintercept=hd05$hd, colour="black", linetype="solid",
               size=2, alpha=0.5) + # thicker median
    geom_vline(data=cdat, aes(xintercept=deciles),
               linetype="solid", size=1, alpha=.5, colour="black") +
    geom_rug() +
    theme_bw() +
    theme(legend.position="none",
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold")) +
    ylab("Density")
  p
}

#' Scatterplots for 2 groups **KEEP?**
#'
#' \code{plot_scat2} produces scatterplots for 2 marginal distributions.
#' The scatterplots are jittered using \code{\link[ggforce]{geom_sina}}.
plot_scat2_sina <- function(data = df,
                       symb_size = 2,
                       symb_stroke = 1,
                       symb_shape = c(21,21),
                       symb_alpha = .2,
                       symb_col = c("black","black"),
                       symb_fil = c("grey70","grey70"),
                       xlabel = NULL,
                       ylabel = NULL,
                       binwidth = NULL,
                       bins = NULL,
                       maxwidth = NULL,
                       scale = TRUE){
  xplot = names(data)[1]
  yplot = names(data)[2]
  p <- ggplot(data, aes_string(x = xplot, y = yplot, fill = xplot,
                               colour = xplot, shape = xplot))
  p <- p + ggforce::geom_sina(size = symb_size, stroke = symb_stroke,
                              alpha = symb_alpha, binwidth = binwidth,
                              bins = bins, maxwidth = maxwidth,
                              scale = scale) +
    theme_bw() +
    scale_colour_manual(values = symb_col) +
    scale_fill_manual(values = symb_fil) +
    scale_shape_manual(values = symb_shape) +
    theme(legend.position="none") +
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
#' @param ... To pass other ggplot2 options.
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

#' Plot deciles and confidence intervals
#'
#' Plot deciles and their confidence intervals
#' using the output of `quantiles_pbci`
#'
#' @export
plot_dec_ci <- function(out = out,
                        plotzero = TRUE,
                        xtitle = "Differences",
                        hjust = -.05,
                        vjust = .2,
                        size_text = 6,
                        colour_dec = "#009E73",
                        fill_dec = "white",
                        size_dec = 4,
                        shape_dec = 21,
                        colour_line = "#009E73",
                        size_line = 1,
                        linetype_line = 2,
                        colour_zero = "black",
                        size_zero = .5,
                        linetype_zero = 1){
  md <- out$est_q[5] # median
  md.c <- as.character(round(md, digits=1)) # turn into characters
  lo.c <- as.character(round(out$ci.low[5], digits=1)) # turn into characters
  up.c <- as.character(round(out$ci.up[5], digits=1)) # turn into characters
  caption <- paste("Median = \n ",md.c," [",lo.c,", ",up.c,"]",sep="")
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
    geom_point(colour = colour_dec,
               size = size_dec,
               shape = shape_dec,
               fill = fill_dec) +
    theme_bw() +
    labs(x="Deciles") +
    labs(y=xtitle) +
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold")) +
    scale_x_continuous(breaks=seq(1,9,1)) +
    annotate("text", x = 5, y = out$ci.up[5], label = caption[1],
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
#' @param xcol,ycol Names of columns of paired observations.
#' @param colour_p Colour parameter of the scatterplot.
#' @param fill_p Fill parameter of the scatterplot and so on.
#' @param colour_q Colour of the segments marking the quartiles.
#' @param linetype_q Linetype of the segments marking the quartiles and so on.
#' @examples
#' plot_scat2d(df, xcol='c1', ycol='c2') # basic call
#' plot_scat2d(df, xcol='c1', ycol='c2', size_q=1) # specify size of quartile segments
#' plot_scat2d(df, xcol='c1', ycol='c2', size_q=c(1,2,1)) # use thicker line for median
#' plot_scat2d(df, xcol='c1', ycol='c2', linetype_q = "longdash") # specify linetype - default = dashed
#' @seealso \code{\link{hd}}
#' @export
plot_scat2d <- function(df = df,
                        xcol = "condition1",
                        ycol = "condition2",
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
  # make data.frames for plotting quartile segments
  hd1.25<-hd(df[[xcol]],.25)
  hd1.5<-hd(df[[xcol]],.5)
  hd1.75<-hd(df[[xcol]],.75)
  hd2.25<-hd(df[[ycol]],.25)
  hd2.5<-hd(df[[ycol]],.5)
  hd2.75<-hd(df[[ycol]],.75)
  df.25<-data.frame(hd1=hd1.25,hd2=hd2.25)
  df.5<-data.frame(hd1=hd1.5,hd2=hd2.5)
  df.75<-data.frame(hd1=hd1.75,hd2=hd2.75)

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
    min.x <- min(df[,2],df[,3])
  }
  if (is.na(max.x)){
    max.x <- max(df[,2],df[,3])
  }
  if (is.na(min.y)){
    min.y <- min(df[,2],df[,3])
  }
  if (is.na(max.y)){
    max.y <- max(df[,2],df[,3])
  }

  # scatterplot of paired observations -----------------
  p <- ggplot(df, aes_string(x=xcol,y=ycol)) +
    # reference line
    geom_abline(intercept = 0) +
    # quartiles
    scale_x_continuous(limits=c(floor(min.x), ceiling(max.x)),breaks=seq(floor(min.x),ceiling(max.x),axis.steps)) +
    scale_y_continuous(limits=c(floor(min.y), ceiling(max.y)),breaks=seq(floor(min.y),ceiling(max.y),axis.steps)) +
    geom_segment(aes(x=hd1,y=min.y,xend=hd1,yend=hd2),data=df.25,linetype=linetype_q[1],size=size_q[1],alpha=alpha_q[1],colour=colour_q[1]) +
    geom_segment(aes(x=hd1,y=min.y,xend=hd1,yend=hd2),data=df.5,linetype=linetype_q[2],size=size_q[2],alpha=alpha_q[2],colour=colour_q[2]) +
    geom_segment(aes(x=hd1,y=min.y,xend=hd1,yend=hd2),data=df.75,linetype=linetype_q[3],size=size_q[3],alpha=alpha_q[3],colour=colour_q[3]) +
    geom_segment(aes(x=min.x,y=hd2,xend=hd1,yend=hd2),data=df.25,linetype=linetype_q[1],size=size_q[1],alpha=alpha_q[1],colour=colour_q[1]) +
    geom_segment(aes(x=min.x,y=hd2,xend=hd1,yend=hd2),data=df.5,linetype=linetype_q[2],size=size_q[2],alpha=alpha_q[2],colour=colour_q[2]) +
    geom_segment(aes(x=min.x,y=hd2,xend=hd1,yend=hd2),data=df.75,linetype=linetype_q[3],size=size_q[3],alpha=alpha_q[3],colour=colour_q[3]) +
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
