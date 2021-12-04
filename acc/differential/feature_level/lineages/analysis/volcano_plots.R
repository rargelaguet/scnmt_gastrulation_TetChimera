library(RColorBrewer)

###########################
## Load default settings ##
###########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera/settings.R")
  source("/Users/ricard/scnmt_gastrulation_TetChimera/acc/differential/feature_level/lineages/analysis/utils.R")
} else {
  stop("Computer not recognised")
}

io$indir <- paste0(io$basedir,"/acc/results/differential/feature_level/lineages")
io$outdir <- paste0(io$basedir,"/acc/results/differential/feature_level/lineages/pdf")

#####################
## Define options ##
#####################

opts$comparisons <- c(
  # "ICM_vs_TE_mural"
  "Prelineage_vs_TE_mural",
  "Prelineage_vs_ICM"
)

# Select genomic contexts
opts$annos <- NULL
if (is.null(opts$annos)) {
  opts$annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
}

# Minimum differential levels (%) for statistical significance
opts$min.diff <- 25

# Multiple testing correction
opts$threshold_fdr <- 0.10

###############
## Load data ##
###############

# Load precomputed differential results
diff.results <- lapply(opts$comparisons, function(i) 
  lapply(opts$annos, function(j) {
    file <- sprintf("%s/%s_%s.txt.gz",io$indir,i,j)
    if (file.exists(file)) {
      fread(file) %>% .[,anno:=as.factor(j)]
    } else {
      cat(sprintf("%s does not exist",file))
    }
  }
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist# %>% .[complete.cases(.)]

# Define statistically significant hits
diff.results[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>opts$min.diff)]

#################
## Filter data ##
#################

# Filter genomic contexs with very few hits
diff.results <- diff.results[,N:=.N,by=c("anno","comparison")] %>% .[N>1000] %>% .[,N:=NULL]


###################
## Volcano plots ##
###################

for (i in unique(diff.results$comparison)) {
  groups <- unlist(strsplit(i,"_vs_"))
  for (j in unique(diff.results$anno)) {
  
    p <- gg_volcano_plot(diff.results[comparison==i & anno==j], groups = groups, top_genes = 0)
    
    # pdf(sprintf("%s/volcano_%s_%s.pdf",io$outdir,i,j), width=7, height=5)
    png(sprintf("%s/volcano_%s_%s.png",io$outdir,i,j), width=7, height=5, units="in", res=400)
    print(p)
    dev.off()
  }
}
