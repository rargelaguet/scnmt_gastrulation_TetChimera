################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/rna/mapping/run/mnn/mapping_mnn.R"
  io$Rscript <- "/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript"
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/rna/mapping/run/mnn/mapping_mnn.R"
  io$tmpdir <- paste0(io$basedir,"/results/mapping/tmp"); dir.create(io$tmpdir)
  stop()
  # io$Rscript <- ""
}
io$outdir <- paste0(io$basedir,"/results/rna/mapping")

####################
## Define options ##
####################

opts$atlas_stages <- c(
  # "E6.5",
  # "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)

opts$query_batches <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

# Test mode (subsetting cells)?
opts$test_mode <- TRUE

#############################
## Run one batch at a time ##
#############################

for (i in opts$query_batches) {
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M 15000 -n 1 -o %s/%s.txt", io$tmpdir,paste(i,collapse=" "))
  }
  cmd <- sprintf("%s %s %s --atlas_stages %s --query_batches %s --outdir %s", lsf, io$Rscript, io$script, paste(opts$atlas_stages,collapse=" "), paste(i,collapse=" "),io$outdir)
  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")
  
  # Run
  print(cmd)
  system(cmd)
}

######################################
## Run all batches at the same time ##
######################################

# Define LSF command
if (grepl("ricard",Sys.info()['nodename'])) {
  lsf <- ""
} else if (grepl("ebi",Sys.info()['nodename'])) {
  lsf <- sprintf("bsub -M 15000 -n 1 -o %s/%s.txt", io$tmpdir,paste(i,collapse=" "))
}
cmd <- sprintf("%s %s %s --atlas_stages %s --query_batches %s --outdir %s", lsf, io$Rscript, io$script, paste(opts$atlas_stages,collapse=" "), paste(opts$query_batches,collapse=" "),io$outdir)
if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")

# Run 
print(cmd)
system(cmd)