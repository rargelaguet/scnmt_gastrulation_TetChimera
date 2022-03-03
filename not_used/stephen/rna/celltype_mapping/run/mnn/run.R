
################
## Define I/O ##
################


source(here::here("settings.R"))
io$script <- here::here("rna/celltype_mapping/run/mnn/mapping_mnn_sjc.R")
io$tmpdir <- paste0(io$basedir,"/results/rna/mapping/tmp"); dir.create(io$tmpdir)

io$outdir <- paste0(io$basedir,"/results/rna/mapping")

####################
## Define options ##
####################

opts$atlas_stages_all <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5",
  "mixed_gastrulation"
)

# opts$query_samples <- c(
#   "E8.5_oct20_plate2",
#   "E8.5_oct20_plate4",
#   "E7.5_tet_crispr_plate5",
#   "E7.5_tet_crispr_plate6",
#   "E7.5_tet_chimera_plate3",
#   "E8.5_oct20_plate1",
#   "E8.5_oct20_plate3",
#   "E8.5_oct20_plate5",
#   "E8.5_oct20_plate6",
#   "E8.5_oct20_plate7",
#   "E8.5_oct20_plate8"      
# )

opts$query_samples <- sample_metadata[, unique(plate)]

# opts$samples <- "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001"

# Test mode (subsetting cells)?
opts$test_mode <- FALSE

if (opts$test_mode) {
  opts$memory <- 10000
} else {
  opts$memory <- 30000
}

opts$npcs <- 50
opts$n_neighbours <- 25

#########
## Run ##
#########

for (i in opts$query_samples) {
  
  # set atlas stages based on query stage
  stages <- sample_metadata[plate == i, unique(stage)] %>%
    gsub("E", "", .) %>%
    as.numeric()
  
  opts$atlas_stages <- paste0("E", c(stages - 0.25, stages, stages + 0.25)) %>% 
    unique() %>%
    .[. %in% opts$atlas_stages_all]
  
  # Define LSF command
  if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M %d -n 1 -o %s/%s.txt", opts$memory, io$tmpdir,paste(i,collapse=" "))
  } else {
    lsf <- ""
  }
  
  cmd <- sprintf(
    "%s Rscript %s --atlas_stages %s --query_samples %s --npcs %d --n_neighbours %d --outdir %s",
    lsf, 
    io$script, 
    paste(opts$atlas_stages,collapse=" "), 
    paste(i,collapse=" "), 
    opts$npcs, 
    opts$n_neighbours, 
    io$outdir
    )
  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")

  # Run
  print(cmd)
  system(cmd)
}
