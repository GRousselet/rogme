plot_diff_asym <- function(data=df){
  # Plot difference asymmetry plot using output from
  # Rand Wilcox's qwmwhd or difQpci functions
  # GAR, University of Glasgow, 2016-07-13
  ylim <- max(max(abs(data$ci.up)),max(abs(data$ci.low)))
  ylim <- c(-ylim,ylim)
  p <- ggplot(data, aes_string(x="quantile", y="SUM")) +
    geom_hline(yintercept=0,linetype=2,alpha=0.5) + # x=0 reference line
    geom_linerange(aes(ymin=ci.low, ymax=ci.up), colour="#009E73") +
    geom_line(colour="#009E73", linetype="solid", size=1.5) +
    geom_point(colour="#009E73", size=4, shape=21, fill="white") + #999999
    xlab("Quantiles") +
    ylab("Quantile sum = q + 1-q") +
    theme_bw() +
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16,face="bold"),
          axis.title.y = element_text(size=16,face="bold")) +
    scale_y_continuous(limits = ylim) +
    scale_x_continuous(breaks = data$quantile)
  p
}
