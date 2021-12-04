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

#############
## Options ##
#############

opts$comparisons <- c(
  "Morula_vs_ICM",
  "Morula_vs_TE",
  "ICM_vs_TE",
  "hESC_vs_ICM",
  "hESC_vs_TE",
  "hESC_vs_Morula"
)

# Select genomic contexts
opts$annos <- NULL
if (is.null(opts$annos)) {
  opts$annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
}

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
      cat(sprintf("%s does not exist\n",file))
    }
  }
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist# %>% .[complete.cases(.)]

# Filter
diff.results <- diff.results[,N:=.N,by=c("anno","comparison")] %>% .[N>1000] %>% .[,N:=NULL]

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
