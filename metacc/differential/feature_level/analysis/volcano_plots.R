library(RColorBrewer)

###########################
## Load default settings ##
###########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/human_embryo_multiomics/met/differential/feature_level/lineages/analysis/load_data.R")
} else {
  stop("Computer not recognised")
}

###################
## Volcano plots ##
###################

for (i in unique(diff.results$comparison)) {
  groups <- unlist(strsplit(i,"_vs_"))
  for (j in unique(diff.results$anno)) {
  
    to.plot <- diff.results[comparison==i & anno==j]
    if (nrow(to.plot)>10) {
      p <- gg_volcano_plot(to.plot, groups = groups, top_genes = 0)
      
      # pdf(sprintf("%s/volcano_%s_%s.pdf",io$outdir,i,j), width=7, height=5)
      png(sprintf("%s/volcano_%s_%s.png",io$outdir,i,j), width=7, height=5, units="in", res=400)
      print(p)
      dev.off()
    }
  }
}
