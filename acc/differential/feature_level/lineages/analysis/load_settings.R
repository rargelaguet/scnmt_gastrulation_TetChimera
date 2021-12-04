library(RColorBrewer)

###########################
## Load default settings ##
###########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera/settings.R")
  source("/Users/ricard/scnmt_gastrulation_TetChimera/acc/differential/feature_level/lineages/analysis/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation_TetChimera/settings.R")
  source("/homes/ricard/scnmt_gastrulation_TetChimera/acc/differential/feature_level/lineages/analysis/utils.R")
} else {
  stop("Computer not recognised")
}

io$indir <- paste0(io$basedir,"/acc/results/differential/feature_level/lineages")
io$outdir <- paste0(io$basedir,"/acc/results/differential/feature_level/lineages/pdf")

#############
## Options ##
#############

opts$lineages <- c(
  "hESC",
  # "Zygote", 
  # "2cell", 
  # "4cell", 
  # "8cell", 
  # "Morula", 
  "ICM", 
  "TE"
)

# Select genomic contexts
opts$annos <- NULL
if (is.null(opts$annos)) {
  opts$annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
}
