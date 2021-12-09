library(RColorBrewer)


###############
## Load data ##
###############

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/human_embryo_multiomics/met/differential/feature_level/lineages/analysis/load_data.R")
} else {
  stop("Computer not recognised")
}

##############
## Barplots ##
##############

to.plot <- diff.results %>%
  # .[,.(number_positive_hits=sum(sig==T & diff>0), number_negative_hits=sum(sig==T & diff<0)), by=c("anno","comparison")] %>%
  .[,.(number_positive_hits=mean(sig==T & diff>0), number_negative_hits=mean(sig==T & diff<0)), by=c("anno","comparison")] %>%
  .[,number_negative_hits:=-number_negative_hits] %>%
  melt(id.vars=c("anno","comparison"))

ylim <- c(min(to.plot$value), max(to.plot$value))

for (i in unique(diff.results$comparison)) {
  p <- gg_barplot(to.plot[comparison==i], title=i, ylim=ylim)
  
  pdf(sprintf("%s/barplots_%s.pdf",io$outdir,i), width=18, height=7)
  # png(sprintf("%s/barplots_%s.png",io$outdir,i), width=6, height=4, units="in", res=400)
  print(p)
  dev.off()
}
