barplot.pub <- function(df, x, colors=NULL, xlim.max=NULL) {
  p <- ggplot(df, aes_string(x=x, y="N")) +
    scale_x_discrete(drop=FALSE) + 
    labs(y="Number of cells") +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.3)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.3), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  if (is.null(colors)) {
    p <- p + geom_bar(stat="identity", color="black")
  } else {
    p <- p + geom_bar(aes_string(fill=x), stat="identity", color="black") + 
      scale_fill_manual(values=colors, drop=F)
  }
  
  if (!is.null(xlim.max)) {
    p <- p + coord_flip(ylim=c(0,xlim.max))
  } else {
    p <- p + coord_flip()
  }
  
  return(p)
}
