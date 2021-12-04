#####################
## Define settings ##
#####################

io <- list(); opts <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$standard.mnn.script <- "/Users/ricard/scnmt_gastrulation_TetChimera/rna/iterative_mapping/run/standard_mnn.R"
  io$iterative.mnn.script <- "/Users/ricard/scnmt_gastrulation_TetChimera/rna/iterative_mapping/run/iterative_mnn.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  io$standard.mnn.script <- "/homes/ricard/scnmt_gastrulation_TetChimera/rna/iterative_mapping/run/standard_mnn.R"
  io$iterative.mnn.script <- "/homes/ricard/scnmt_gastrulation_TetChimera/rna/iterative_mapping/run/iterative_mnn.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_gastrulation_TetChimera/rna/results/iterative_mapping/tmp"
} 

# Atlas stages
opts$atlas_stages <- c(
  # "E6.5"
  # "E6.75",
  # "E7.0",
  "E7.25",
  "E7.5",
  "E7.75"
  # "E8.0",
  # "E8.25",
  # "E8.5"
  # "mixed_gastrulation"
)

# Query plates
# opts$query_plates <- c(
#   # "tet_chimera_march20_plate1",
#   # "tet_chimera_march20_plate2",
#   # "tet_chimera_march20_plate3",
#   # "tet_chimera_march20_plate4"
#   "tet_chimera_march20_plate5",
#   "tet_chimera_march20_plate6"
# )

opts$query_plates <- c(
  "tet_chimera_oct19_plate1",
  "tet_chimera_oct19_plate2",
  "tet_chimera_oct19_plate4",
  "tet_chimera_oct19_plate5"
)

# Test mode (subset cells)?
opts$test <- FALSE


#################################
## Run each plate individually ##
#################################

for (i in opts$query_plates) {
  # LSF
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else {
    lsf <- sprintf("bsub -M 30000 -n 1 -q research-rh74 -o %s/%s.txt", io$tmpdir, i)
  }

  # Run standard MNN
  cmd <- sprintf("%s Rscript %s --query_plates %s --atlas_stages %s", lsf, io$standard.mnn.script, i, paste(opts$atlas_stages, collapse=" "))
  if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
  system(cmd)

  # Run tree-guided MNN
  cmd <- sprintf("%s Rscript %s --query_plates %s --atlas_stages %s", lsf, io$iterative.mnn.script, i, paste(opts$atlas_stages, collapse=" "))
  if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
  system(cmd)
}


#############################
## Run all plates together ##
#############################

# LSF
if (grepl("ricard",Sys.info()['nodename'])) {
  lsf <- ""
} else {
  lsf <- sprintf("bsub -M 30000 -n 1 -q research-rh74 -o %s/all_plates.txt", io$tmpdir)
}

# Run standard MNN
cmd <- sprintf("%s Rscript %s --query_plates %s --atlas_stages %s", lsf, io$standard.mnn.script, paste(opts$query_plates, collapse=" "), paste(opts$atlas_stages, collapse=" "))
if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
system(cmd)

# Run tree-guided MNN
cmd <- sprintf("%s Rscript %s --query_plates %s --atlas_stages %s", lsf, io$iterative.mnn.script, paste(opts$query_plates, collapse=" "), paste(opts$atlas_stages, collapse=" "))
if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
system(cmd)
