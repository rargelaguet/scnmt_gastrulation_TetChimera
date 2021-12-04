###########################
## Load default settings ##
###########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera//settings.R")
  io$script <- "/Users/ricard/scnmt_gastrulation_TetChimera/met/differential/feature_level/pseudobulk/diffmet_feature_level.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/scnmt_gastrulation_TetChimera//settings.R")
  # io$script <- "/homes/ricard/scnmt_gastrulation_TetChimera/met/differential/feature_level/diffmet_feature_level.R"
} else {
  stop("Computer not recognised")
}

io$outdir <- paste0(io$basedir,"/met/differential/feature_level/pseudobulk")

#############
## Options ##
#############

# opts$groups <- list(
#   "ICM_vs_TE_mural" = list(c("ICM"), c("TE_mural")),
#   "ICM_vs_TE_polar" = list(c("ICM"), c("TE_polar"))
# )


# Genomic contexts
opts$annos <- NULL
if (is.null(opts$annos)) {
  opts$annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
}
# opts$annos <- c("prom_2000_2000")

#########
## Run ##
#########

# for (i in names(opts$groups)) {
for (i in 1:(length(opts$lineages)-1)) {
  groupA <- opts$lineages[[i]]
  for (j in (i+1):length(opts$lineages)) {
    groupB <- opts$lineages[[j]]
    for (anno in opts$annos) {
      outfile <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$outdir, groupA, groupB, anno)
      # lsf <- sprintf("bsub -M 8000 -n 1 -q research-rh74 -o %s/%s_%s.txt", io$tmpdir, i, j)
      lsf <- ""
      cmd <- sprintf("%s Rscript %s --anno %s --groupA %s --groupB %s --outfile %s", 
                     lsf, io$script, anno, groupA, groupB, outfile)
      print(cmd)
      system(cmd)
    }
  }
}
