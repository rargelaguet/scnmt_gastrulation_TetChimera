here::here("metacc/mofa/load_settings.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outfile <- paste0(io$basedir,"/results/metacc/mofa/model.rds"); dir.create(dirname(io$outfile), showWarnings = F)

# Reticulate
io$python <- "/Users/rargelaguet/opt/anaconda3/envs/main/bin/python"
reticulate::use_python(io$python, required = T)

# Define genomic annotations
opts$met_annos <- c("multiome_peaks")
opts$acc_annos <- c("multiome_peaks")

# Define which stage and lineages to look at 
# opts$lineages <- c(
#   "E5.5_Epiblast",
#   "E5.5_Visceral_endoderm"
# )

# Filtering options for methylation
opts$met_min_observations <- 1        # minimum number of CpG sites per feature
opts$met_min_cells <- 50      # minimum number of cells per feature (per stage)
opts$met_nfeatures <- 2500    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min_observations <- 5        # minimum number of GpC sites per feature
opts$acc_min_cells <- 50      # minimum number of cells per feature (per stage)
opts$acc_nfeatures <- 2500    # maximum number of features per view (filter based on variance)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[pass_metQC==TRUE & pass_accQC==TRUE] %>%
  setnames("celltype.mapped","celltype")

opts$met_cells <- cell_metadata.dt %>% .[pass_metQC==T,id_met]
opts$acc_cells <- cell_metadata.dt %>% .[pass_accQC==T,id_acc]
