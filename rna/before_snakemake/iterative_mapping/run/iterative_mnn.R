# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(edgeR))

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

########################################################
## Define distance matrix for hierarchical clustering ##
########################################################

opts$celltypes <- which(table(meta_atlas$celltype)>25) %>% names# %>% stringr::str_replace_all("_", " ")
# opts$celltypes <- unique(sample_metadata_atlas$celltype) 

dist <- fread(paste0(io$atlas.basedir,"/results/phylogenetic_tree/PAGA_distances.csv.gz")) %>%
  as.data.frame %>% tibble::column_to_rownames("V1") %>% as.matrix %>%
  .[opts$celltypes,opts$celltypes] %>%
  as.dist

#######################
## Recursive mapping ##
#######################

sce.query$celltype_mapped <- paste(opts$celltypes,collapse="%")

while (any(grepl("%",sce.query$celltype_mapped))) {
  print(table(sce.query$celltype_mapped))
  mapping_dt <- recursive.fn(sce.query, sce.atlas, dist)
  ids <- match(mapping_dt$cell,colnames(sce.query))
  sce.query$celltype_mapped[ids] <- mapping_dt$celltype_mapped
  sce.query$mapping_score[ids] <- mapping_dt$mapping_score
}

##########
## Save ##
##########

mapping_dt <- data.table(
  cell = colnames(sce.query), 
  celltype_mapped = sce.query$celltype_mapped,
  mapping_score = sce.query$mapping_score
)
# foo <- mapping_dt %>% merge(meta_query[,c("cell","celltype.mapped","celltype.score")] %>% setnames(c("cell","celltype_old","score_old")))

fwrite(mapping_dt, sprintf("%s/%s_iterative_mnn.txt.gz",io$outdir,paste(args$query_plates,collapse="_")), sep="\t")

