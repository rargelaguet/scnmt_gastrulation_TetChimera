#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera//settings.R")
  io$script <- "/Users/ricard/scnmt_gastrulation_TetChimera//met/differential/feature_level/lineages/diffmet_feature_level.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/scnmt_gastrulation_TetChimera//settings.R")
  io$script <- "/homes/ricard/scnmt_gastrulation_TetChimera//met/differential/feature_level/lineages/diffmet_feature_level.R"
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/met/results/differential/feature_level/lineages"); dir.create(io$outdir, showWarnings=F, recursive = T)
io$tmpdir <- paste0(io$basedir,"/met/results/differential/feature_level/lineages/tmp")

#############
## Options ##
#############

# Minimum number of cells per group
opts$min.cells <- 10

# Define lineages
opts$lineages <- c(
  "hESC",
  "Zygote", 
  "2cell",
  "4cell", 
  "8cell", 
  "Morula", 
  "ICM", 
  "TE"
)

# Genomic contexts
# opts$annos <- c("prom_2000_2000")
# opts$annos <- NULL
if (is.null(opts$annos)) {
  opts$annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
  # io$diff.met <- "/Users/ricard/data/scnmt_gastrulation_TetChimera/met/feature_level"
  # opts$annos <- list.files(io$diff.met, pattern=".tsv.gz") %>% gsub(".tsv.gz","",.)
}

#########
## Run ##
#########

for (i in 1:(length(opts$lineages)-1)) {
  groupA <- opts$lineages[[i]]
  for (j in (i+1):length(opts$lineages)) {
    groupB <- opts$lineages[[j]]
    for (anno in opts$annos) {
      outfile <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$outdir, groupA, groupB, anno)
      if (file.exists(outfile)) {
        print(sprintf("%s already exists, skipping...",outfile))
      } else {
        lsf <- sprintf("bsub -M 8000 -n 1 -q research-rh74 -o %s/%s_%s_%s.txt", io$tmpdir, groupA, groupB, anno)
        # lsf <- ""
        cmd <- sprintf("%s Rscript %s --anno %s --groupA %s --groupB %s --min.cells %d --outfile %s", 
                       lsf, io$script, anno, paste(groupA, collapse=" "), paste(groupB, collapse=" "), opts$min.cells, outfile)
        print(cmd)
        system(cmd)
      }
    }
  }
}


##############################
## Run selected comparisons ##
##############################

# groupA <- c("Prelineage")
# groupB <- c("Morula")

# # for (anno in opts$annos) {
#   anno <- "LINE"
#   outfile <- sprintf("%s/%s_vs_%s_%s.txt.gz", io$outdir, groupA, groupB, anno)
#   lsf <- ""
#   cmd <- sprintf("%s Rscript %s --anno %s --groupA %s --groupB %s --min.cells %d --outfile %s",
#                  lsf, io$script, anno, paste(groupA, collapse=" "), paste(groupB, collapse=" "), opts$min.cells, outfile)
#   print(cmd)
#   system(cmd)
# # }
