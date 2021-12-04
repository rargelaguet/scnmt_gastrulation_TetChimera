#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera//settings.R")
  io$script <- "/Users/ricard/scnmt_gastrulation_TetChimera//met/differential/feature_level/sex/diffmet_feature_level.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/scnmt_gastrulation_TetChimera//settings.R")
  io$script <- "/homes/ricard/scnmt_gastrulation_TetChimera//met/differential/feature_level/sex/diffmet_feature_level.R"
} else {
  stop("Computer not recognised")
}

io$outdir <- paste0(io$basedir,"/met/results/differential/feature_level/sex")
io$tmpdir <- paste0(io$basedir,"/met/results/differential/feature_level/sex/tmp")

# table(sample_metadata[pass_metQC==T,lineage])

#############
## Options ##
#############

# Minimum number of cells per group
opts$min.cells <- 10

# Define lineages
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

# Genomic contexts
# opts$annos <- c("genebody")
opts$annos <- NULL
if (is.null(opts$annos)) {
  # opts$annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
  # io$diff.met <- "/Users/ricard/data/scnmt_gastrulation_TetChimera/met/feature_level"
  opts$annos <- list.files(io$met_data_parsed, pattern=".tsv.gz") %>% gsub(".tsv.gz","",.)
}

#########
## Run ##
#########

for (i in opts$lineages) {
  for (j in opts$annos) {
    outfile <- sprintf("%s/%s_%s_Male_vs_Female.txt.gz", io$outdir,j,i)
    if (file.exists(outfile)) {
      print(sprintf("%s already exists, skipping...",outfile))
    } else {
      lsf <- sprintf("bsub -M 15000 -n 1 -q research-rh74 -o %s/%s_%s_Male_vs_Female.txt", io$tmpdir, i, j)
      # lsf <- ""
      cmd <- sprintf("%s Rscript %s --anno %s --lineage %s --min.cells %d --outfile %s", lsf, io$script, j, i, opts$min.cells, outfile)
      print(cmd)
      system(cmd)
    }
  }
}
