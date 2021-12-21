#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
  io$script <- "/Users/ricard/10x_gastrulation_TetChimera/differential/differential.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")
  io$script <- "/homes/ricard/10x_gastrulation_TetChimera/differential/differential.R"
  io$tmpdir <- paste0(io$basedir,"/results/differential/tmp"); dir.create(io$tmpdir, showWarnings = F)
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/differential"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Statistical test
opts$statistical.test <- "edgeR"

# Testing mode
opts$test_mode <- FALSE

# Define groups
opts$groupA <- c("E8.5_Host")
opts$groupB <- c("E8.5_TET_TKO")

#########
## Run ##
#########

# for (i in head(opts$celltypes,n=3)) {
for (i in opts$celltypes) {
    outfile <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$outdir,i,opts$groupA,opts$groupB)
    
    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M 15000 -n 1 -o %s/%s_%s_vs_%s.txt", io$tmpdir,i,opts$groupA,opts$groupB)
    }
    cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --celltype %s --test %s --outfile %s", lsf, io$script, opts$groupA, opts$groupB, i, opts$statistical.test, outfile)
    if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
    
    # Run
    print(cmd)
    system(cmd)
}

