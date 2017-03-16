#' Plot 2 groups: KDE + rug + deciles
#'
#' Plot kernel density estimates + rug plots + deciles
#' for 2 groups stored in a data frame.
#'
#' @export
plot_kde_rug_dec2 <- function(data = df){
  cdat <- plyr::ddply(data, "gr", summarise, deciles = q1469(data))
  hd05 <- plyr::ddply(data, "gr", summarise, hd = hd(data,0.5))
  #cc <- "grey80" # colour to plot deciles
  p <- ggplot(data, aes(x=data, fill=gr)) + geom_density(alpha=.3) +
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
  cdat <- plyr::ddply(data, "gr", summarise, deciles=q1469(data))
  hd05 <- plyr::ddply(data, "gr", summarise, hd=hd(data,0.5))
  cc <- "grey80" # colour to plot deciles
  p <- ggplot(data, aes(x=data)) +
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

#' Scatterplots for 2 groups
#'
#' \code{plot_scat2} produces scatterplots for 2 marginal distributions. #' The scatterplots are jittered using \code{\link[ggforce]{geom_sina}}.
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

#' Scatterplots for 2 groups
#'
#' \code{plot_scat2} produces scatterplots for 2 marginal distributions.
#' The scatterplots are jittered using \code{\link[ggbeeswarm]{geom_quasirandom}}.
#' @export
plot_scat2 <- function(data = df,
                       xlabel = NULL,
                       ylabel = NULL,
                       ...){
  xplot = names(data)[1]
  yplot = names(data)[2]
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

#' Plot deciles and their confidence intervals
#' using the output of `quantiles_pbci`
#' GAR, University of Glasgow, 2016-07-15
#' @export
plot_dec_ci <- function(out = out,
                        plotzero = TRUE,
                        xtitle = "Differences",
                        hjust = -.05,
                        vjust = .2,
                        size_text = 6,
                        colour_dec = "#009E73",
                        fill_dec = "white",
                        colour_line = "#009E73"){
  md <- out$est_q[5] # median
  md.c <- as.character(round(md, digits=1)) # turn into characters
  lo.c <- as.character(round(out$ci.low[5], digits=1)) # turn into characters
  up.c <- as.character(round(out$ci.up[5], digits=1)) # turn into characters
  caption <- paste("Median = \n ",md.c," [",lo.c,", ",up.c,"]",sep="")
  p <- ggplot(data=out, aes(x=quantile*10, y=est_q)) +
    geom_abline(intercept = md, slope = 0,colour="black",size=.5,linetype=2)
  if (plotzero){
    p <- p + geom_abline(intercept = 0, slope = 0,colour="black",size=.5,linetype=1)
  }
  p <- p + geom_linerange(aes(ymin=ci.low, ymax=ci.up), colour=colour_line,size=1) +
    geom_point(colour=colour_dec, size=4, shape=21, fill=fill_dec) +
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
#' Quartiles of each condition are superimposed.
#' Input is a data frame with 3 columns: participant, condition1, condition 2
#' @export
plot_scat2d <- function(df=df,
                        xname="condition1",
                        yname="condition2",
                        min.x=NA,
                        min.y=NA,
                        max.x=NA,
                        max.y=NA,
                        axis.steps=2,
                        psize=5,
                        pstroke=1,
                        pshape=21,
                        pcolour="black",
                        pfill="#ffb347",
                        palpha=.5){
  # make data.frames for plotting quartile segments
  hd1.5<-hd(df[,2],.5)
  hd1.25<-hd(df[,2],.25)
  hd1.75<-hd(df[,2],.75)
  hd2.5<-hd(df[,3],.5)
  hd2.25<-hd(df[,3],.25)
  hd2.75<-hd(df[,3],.75)
  df.5<-data.frame(hd1=hd1.5,hd2=hd2.5)
  df.25<-data.frame(hd1=hd1.25,hd2=hd2.25)
  df.75<-data.frame(hd1=hd1.75,hd2=hd2.75)

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
  p <- ggplot(df, aes_string(x=xname,y=yname)) +
    geom_abline(intercept = 0) +
    geom_point(size=psize,
               stroke=pstroke,
               shape=pshape,
               colour=pcolour,
               fill=pfill,
               alpha=palpha) +
    theme_bw() +
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold"),
          legend.title = element_blank(),
          plot.title = element_text(size=20)) +
    labs(title="Paired observations") +
    scale_x_continuous(limits=c(floor(min.x), ceiling(max.x)),breaks=seq(floor(min.x),ceiling(max.x),axis.steps)) +
    scale_y_continuous(limits=c(floor(min.y), ceiling(max.y)),breaks=seq(floor(min.y),ceiling(max.y),axis.steps)) +
    geom_segment(aes(x=hd1,y=min.y,xend=hd1,yend=hd2),data=df.5,linetype="dashed",size=1,alpha=.5,colour="black") +
    geom_segment(aes(x=hd1,y=min.y,xend=hd1,yend=hd2),data=df.25,linetype="dashed",size=.5,alpha=.5,colour="black") +
    geom_segment(aes(x=hd1,y=min.y,xend=hd1,yend=hd2),data=df.75,linetype="dashed",size=.5,alpha=.5,colour="black") +
    geom_segment(aes(x=min.x,y=hd2,xend=hd1,yend=hd2),data=df.5,linetype="dashed",size=1,alpha=.5,colour="black") +
    geom_segment(aes(x=min.x,y=hd2,xend=hd1,yend=hd2),data=df.25,linetype="dashed",size=.5,alpha=.5,colour="black") +
    geom_segment(aes(x=min.x,y=hd2,xend=hd1,yend=hd2),data=df.75,linetype="dashed",size=.5,alpha=.5,colour="black")

  #   geom_vline(xintercept=hd(df$condition1,.5),linetype="dashed", size=1, alpha=.3, colour="black") +
  #   geom_vline(xintercept=hd(df$condition1,.25),linetype="dashed", size=.5, alpha=.3, colour="black") +
  #   geom_vline(xintercept=hd(df$condition1,.75),linetype="dashed", size=.5, alpha=.3, colour="black") +
  #   geom_hline(yintercept=hd(df$condition2,.5),linetype="dashed", size=1, alpha=.3, colour="black") +
  #   geom_hline(yintercept=hd(df$condition2,.25),linetype="dashed", size=.5, alpha=.3, colour="black") +
  #   geom_hline(yintercept=hd(df$condition2,.75),linetype="dashed", size=.5, alpha=.3, colour="black")
  p
}
# ----------------------------------------------------------------------------
