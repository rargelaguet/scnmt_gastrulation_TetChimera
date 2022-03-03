suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))

source(here::here("settings.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_samples',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--query_sce',       type="character",               help='SingleCellExperiment file for the query')
p$add_argument('--atlas_sce',       type="character",               help='SingleCellExperiment file for the atlas')
p$add_argument('--query_metadata',  type="character",               help='metadata file for the query')
p$add_argument('--atlas_metadata',  type="character",               help='metadata file for the atlas')
p$add_argument('--npcs',            type="integer",                 help='Number of principal components')
p$add_argument('--n_neighbours',    type="integer",                 help='Number of neighbours')
p$add_argument('--test',            action = "store_true",          help='Testing mode')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))


args$query_sce <- io$sce
args$query_metadata <- io$metadata
args$atlas_metadata <- io$atlas.metadata
args$atlas_sce <- io$atlas.sce


#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna/mapping/run/mapping_functions.R")
} else {
  source(here::here("settings.R"))
  source(here::here("rna/celltype_mapping/run/mapping_functions.R"))
}
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir

## START TEST ##

# args$atlas_metadata <- io$atlas.metadata
# args$atlas_sce <- io$atlas.sce
# opts$atlas_stages <- c(
#   # "E6.5",
#   # "E6.75",
#   # "E7.0",
#   "E7.25",
#   "E7.5",
#   "E7.75"
#   "E8.0",
#   "E8.25",
#   "E8.5"
#   # "mixed_gastrulation"
# )
# args$n_neighbours <- 25
# args$npcs <- 50
# args$outdir <- paste0(io$basedir, "/rna/results/mapping")
# args$query_metadata <- io$metadata
# args$query_samples <- c(
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
# args$query_sce <- io$sce
# args$test <- FALSE


## END TEST ##

if (isTRUE(args$test)) print("Test mode activated...")

####################
## Define options ##
####################

################
## Load atlas ##
################

# Load cell metadata
meta_atlas <- fread(args$atlas_metadata) %>%
  .[stripped==F & doublet==F & stage%in%args$atlas_stages]

# Filter
if (isTRUE(args$test)) meta_atlas <- head(meta_atlas,n=1000)

# Load SingleCellExperiment
sce_atlas <- load_SingleCellExperiment(
  args$atlas_sce, normalise = TRUE, 
  cells = meta_atlas$cell, 
  remove_non_expressed_genes = TRUE
  )

if (!"cell"%in%colnames(colData(sce_atlas))) {
  sce_atlas$cell <- colnames(sce_atlas)
}

################
## Load query ##
################

# Load cell metadata
meta_query <- fread(args$query_metadata) %>% 
  .[, c("sample", "cell") := .(plate, id_rna)] %>%        # reformat names to match 10x style
  .[pass_rnaQC == TRUE & sample %in% args$query_samples] 
                                        
if (isTRUE(args$test)) meta_query <- head(meta_query,n=1000)

# Load SingleCellExperiment
sce_query <- load_SingleCellExperiment(args$query_sce, cells = meta_query$cell, remove_non_expressed_genes = TRUE)

if (!"cell"%in%colnames(colData(sce_query))) {
  sce_query$cell <- colnames(sce_query)
}

#############
## Prepare ## 
#############

# Rename ensemble IDs to gene names in the atlas
gene_metadata <- fread(io$gene_metadata) %>% .[,c("chr","ens_id","symbol")] %>%
  .[symbol!="" & ens_id%in%rownames(sce_atlas)] %>%
  .[!duplicated(symbol)]

sce_atlas <- sce_atlas[rownames(sce_atlas)%in%gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
rownames(sce_atlas) <- foo[rownames(sce_atlas)]

# Sanity cehcks
stopifnot(sum(is.na(rownames(sce_atlas)))==0)
stopifnot(sum(duplicated(rownames(sce_atlas)))==0)

# Rename ensemble IDs to gene names in query
gene_metadata <- fread(io$gene_metadata) %>% .[,c("chr","ens_id","symbol")] %>%
  .[symbol!="" & ens_id%in%rownames(sce_query)] %>%
  .[!duplicated(symbol)]

sce_query <- sce_query[rownames(sce_query)%in%gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
rownames(sce_query) <- foo[rownames(sce_query)]

# Sanity cehcks
stopifnot(sum(is.na(rownames(sce_query)))==0)
stopifnot(sum(duplicated(rownames(sce_query)))==0)

# Intersect genes
genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))

# Filter some genes manually
genes.intersect <- genes.intersect[grep("^mt-",genes.intersect,invert = T)]
genes.intersect <- genes.intersect[grep("^Rik",genes.intersect,invert = T)]
genes.intersect <- genes.intersect[grep("Rps|Rpl",genes.intersect,invert = T)]
genes.intersect <- genes.intersect[!genes.intersect=="Xist"]
genes.intersect <- genes.intersect[!genes.intersect%in%gene_metadata[chr=="chrY",symbol]]

# Subset SingleCellExperiment objects
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

#########
## Map ##
#########

mapping  <- mapWrap(
  atlas_sce = sce_atlas, 
  atlas_meta = meta_atlas,
  map_sce = sce_query, 
  map_meta = meta_query, 
  # genes = marker_genes, 
  npcs = args$npcs, 
  k = args$n_neighbours
)

##########
## Save ##
##########

outfile <- sprintf("%s/mapping_mnn_%s.rds",args$outdir,paste(args$query_samples,collapse="-"))
dir.create(args$outdir, recursive = TRUE)
#outfile <- paste0(args$outdir, "/mapping_mnn.rds")
print(paste("saving as:", outfile))
saveRDS(mapping, outfile)

##########
## TEST ##
##########

# foo <- mapping$mapping %>% as.data.table() %>% .[,c("cell","celltype.mapped","celltype.score")]
# bar <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/rna/mapping/sample_metadata_after_mapping.txt.gz") %>% 
#   .[,c("cell","celltype.mapped","celltype.score")]
# foobar <- merge(foo,bar,by=c("cell"))
# foobar[celltype.mapped.x!=celltype.mapped.y] %>% View
# mean(foobar$celltype.mapped.x==foobar$celltype.mapped.y)
# fwrite(foo,"/Users/ricard/data/gastrulation_multiome_10x/results/rna/mapping/E7.5_rep2_test.txt.gz")
