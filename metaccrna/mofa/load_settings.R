here::here("metaccrna/mofa/load_settings.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Python binary for reticulate connection
io$outfile <- paste0(io$basedir,"/results/metaccrna/mofa2/model.rds"); dir.create(dirname(io$outfile), showWarnings = F)

# io$python <- "..."

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
opts$met_min_cells <- 10      # minimum number of cells per feature (per stage)
opts$met_nfeatures <- 2500    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min_observations <- 5        # minimum number of GpC sites per feature
opts$acc_min_cells <- 10      # minimum number of cells per feature (per stage)
opts$acc_nfeatures <- 2500    # maximum number of features per view (filter based on variance)

# Filtering options for RNA
opts$rna_min_cdr <- 0.25      # Remove genes with small cellular detection rate
opts$rna_ngenes <- 2500       # maximum number of genes (filter based on variance)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE | pass_metQC==TRUE | pass_accQC==TRUE] %>%
  setnames("celltype.mapped","celltype")

opts$met_cells <- cell_metadata.dt %>% .[pass_metQC==T,id_met]
opts$acc_cells <- cell_metadata.dt %>% .[pass_accQC==T,id_acc]
opts$rna_cells <- cell_metadata.dt %>% .[pass_rnaQC==T,id_rna]

