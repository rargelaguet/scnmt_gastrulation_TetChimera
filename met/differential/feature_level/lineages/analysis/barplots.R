library(RColorBrewer)

###########################
## Load default settings ##
###########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera/public_data/Li2016/settings.R")
  source("/Users/ricard/scnmt_gastrulation_TetChimera/public_data/Li2016/met/differential/feature_level/analysis/utils.R")
} else {
  stop("Computer not recognised")
}

io$indir <- paste0(io$basedir,"/met/differential/feature_level")
io$outdir <- paste0(io$basedir,"/met/differential/feature_level/pdf")

#############
## Options ##
#############

opts$comparisons <- c(
  "ICM_vs_TE"
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
      cat(sprintf("%s does not exist",file))
    }
  }
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist# %>% .[complete.cases(.)]

# Filter
diff.results <- diff.results[,N:=.N,by=c("anno","comparison")] %>% .[N>1000] %>% .[,N:=NULL]

##############
## Barplots ##
##############

tmp <- diff.results %>%
  .[,.(number_positive_hits=mean(sig==TRUE & diff>0), number_negative_hits=mean(sig==TRUE & diff<0)), by=c("anno","comparison")] %>%
  .[,number_negative_hits:=-number_negative_hits] %>%
  melt(id.vars=c("anno","comparison"))

ylim <- c(min(tmp$value), max(tmp$value))

for (i in unique(diff.results$comparison)) {
  p <- gg_barplot(tmp[comparison==i], title=i, ylim=ylim)
  
  pdf(sprintf("%s/barplots_%s.pdf",io$outdir,i), width=18, height=7)
  # png(sprintf("%s/barplots_%s.png",io$outdir,i), width=6, height=4, units="in", res=400)
  print(p)
  dev.off()
}
