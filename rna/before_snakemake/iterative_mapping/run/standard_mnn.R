# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_plates',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--test',            action = "store_true",  help = 'Testing mode')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$atlas_stages <- c(
#   # "E6.5"
#   # "E6.75",
#   # "E7.0",
#   "E7.25",
#   "E7.5",
#   "E7.75"
#   # "E8.0",
#   # "E8.25",
#   # "E8.5"
#   # "mixed_gastrulation"
# )
# 
# args$query_plates <- c(
#   # "tet_chimera_march20_plate1"
#   # "tet_chimera_march20_plate2",
#   # "tet_chimera_march20_plate3",
#   # "tet_chimera_march20_plate4",
#   # "tet_chimera_march20_plate5",
#   # "tet_chimera_march20_plate6"
#   "tet_chimera_oct19_plate1",
#   "tet_chimera_oct19_plate2",
#   "tet_chimera_oct19_plate4",
#   "tet_chimera_oct19_plate5"
# )
# 
# args$test <- FALSE
## END TEST ##

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera/settings.R")
  source("/Users/ricard/scnmt_gastrulation_TetChimera/rna/iterative_mapping/run/utils.R")
  io$atlas.marker_genes <- "/Users/ricard/data/gastrulation10x/results/marker_genes/E7.0_to_E7.75/marker_genes.txt.gz"
  io$script_load_data <- "/Users/ricard/scnmt_gastrulation_TetChimera/rna/iterative_mapping/run/load_data.R"
  
} else {
  source("/homes/ricard/scnmt_gastrulation_TetChimera/settings.R")
  source("/homes/ricard/scnmt_gastrulation_TetChimera/rna/iterative_mapping/run/utils.R")
  io$atlas.marker_genes <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/marker_genes/E7.0_to_E7.75/marker_genes.txt.gz"
  io$script_load_data <- "/homes/ricard/scnmt_gastrulation_TetChimera/rna/iterative_mapping/run/load_data.R"
}
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir
io$outdir <- paste0(io$basedir,"/rna/results/iterative_mapping")

if (isTRUE(args$test)) print("Test mode activated...")

###############
## Load data ##
###############

io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_gastrulation_TetChimera/TO_MERGE_scnmt_gastrulation_TetKO"
io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
io$sce <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")

source(io$script_load_data)

#############
## Run MNN ##
#############

mapping_dt <- mnn.fn(sce_query, sce_atlas, npcs = 50, k = 25, cosineNorm = TRUE)

# table(mapping_dt$celltype_mapped)

# foo <- fread("/Users/ricard/data/scnmt_gastrulation_TetChimera/rna/results/iterative_mapping/tet_chimera_march20_plate1_standard_mnn.txt.gz") %>%
#   merge(mapping_dt, by="cell")

# foo <- fread("/Users/ricard/data/scnmt_gastrulation_TetChimera/TO_MERGE_scnmt_gastrulation_TetKO/rna/mapping_10x/sample_metadata_mapping_mnn.txt") %>%
#   merge(mapping_dt %>% copy %>% setnames("cell","id_rna"), by="id_rna") %>%
#   .[,c("id_rna","lineage10x","celltype_mapped")]

##########
## Save ##
##########

fwrite(mapping_dt, sprintf("%s/%s_standard_mnn.txt.gz",io$outdir,paste(args$query_plates,collapse="_")), sep="\t")

