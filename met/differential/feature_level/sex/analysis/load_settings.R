library(RColorBrewer)

###########################
## Load default settings ##
###########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera//settings.R")
  source("/Users/ricard/scnmt_gastrulation_TetChimera//met/differential/feature_level/sex/analysis/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation_TetChimera//settings.R")
  source("/homes/ricard/scnmt_gastrulation_TetChimera//met/differential/feature_level/sex/analysis/utils.R")
} else {
  stop("Computer not recognised")
}

io$indir <- paste0(io$basedir,"/met/results/differential/feature_level/sex")
io$outdir <- paste0(io$basedir,"/met/results/differential/feature_level/sex/pdf")

#############
## Options ##
#############

opts$lineages <- c(
  # "hESC",
  # "Zygote", 
  # "2cell", 
  # "4cell", 
  "8cell", 
  # "Morula", 
  "ICM", 
  "TE"
)

# Select genomic contexts
opts$annos <- NULL
if (is.null(opts$annos)) {
  opts$annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
}
