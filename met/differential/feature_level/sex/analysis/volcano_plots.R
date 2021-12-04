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

###################
## Volcano plots ##
###################

groups <- c("Male","Female")

for (i in unique(diff.results$lineage)) {
  for (j in unique(diff.results$anno)) {
  
    to.plot <- diff.results[lineage==i & anno==j]
    
    p <- gg_volcano_plot(to.plot, groups, top_genes = 0)
    
    # pdf(sprintf("%s/volcano_%s_%s.pdf",io$outdir,i,j), width=7, height=5)
    png(sprintf("%s/volcano_%s_%s.png",io$outdir,i,j), width=7, height=5, units="in", res=400)
    print(p)
    dev.off()
  }
}
