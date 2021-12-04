###########################
## Load default settings ##
###########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera//met/differential/feature_level/sex/analysis/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation_TetChimera//met/differential/feature_level/sex/analysis/load_settings.R")
} else {
  stop("Computer not recognised")
}

###############
## Load data ##
###############

# Load precomputed differential results
diff.results <- lapply(opts$lineages, function(i) 
  lapply(opts$annos, function(j) {
    file <- sprintf("%s/%s_%s_Male_vs_Female.txt.gz",io$indir,j,i)
    if (file.exists(file)) {
      fread(file) %>% .[,anno:=as.factor(j)]
    } else {
      cat(sprintf("%s does not exist",file))
    }
  }
  ) %>% rbindlist %>% .[,lineage:=i] 
) %>% rbindlist

# Filter
diff.results <- diff.results[,N:=.N,by=c("anno","lineage")] %>% .[N>500] %>% .[,N:=NULL]

##############
## Barplots ##
##############

to.plot <- diff.results %>%
  # .[,.(number_positive_hits=sum(sig==T & diff>0), number_negative_hits=sum(sig==T & diff<0)), by=c("anno","lineage")] %>%
  .[,.(number_positive_hits=mean(sig==T & diff>0), number_negative_hits=mean(sig==T & diff<0)), by=c("anno","lineage")] %>%
  .[,number_negative_hits:=-number_negative_hits] %>%
  melt(id.vars=c("anno","lineage"))

ylim <- c(min(to.plot$value), max(to.plot$value))

for (i in unique(diff.results$lineage)) {
  p <- gg_barplot(to.plot[lineage==i], title=i, ylim=ylim)
  
  pdf(sprintf("%s/barplots_%s.pdf",io$outdir,i), width=18, height=7)
  # png(sprintf("%s/barplots_%s.png",io$outdir,i), width=6, height=4, units="in", res=400)
  print(p)
  dev.off()
}
