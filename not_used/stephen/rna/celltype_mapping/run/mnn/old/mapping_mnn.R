# suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_samples',   type="character",   nargs='+',  help='Query sample(es)')
p$add_argument('--test',            action = "store_true",          help='Testing mode')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

# # START TEST ##
args$atlas_stages <- c(
  # "E6.5",
  # "E6.75",
  # "E7.0",
  # "E7.25",
  "E7.5"
  # "E7.75",
  # "E8.0",
  # "E8.25",
  # "E8.5"
  # "mixed_gastrulation"
)
args$query_samples <- c("E7.5_rep1")
args$test <- TRUE
args$outdir <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/mapping"
## END TEST ##

################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna/mapping/run/mapping_functions.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/rna/mapping/run/mapping_functions.R")
}
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir

if (isTRUE(args$test)) print("Test mode activated...")

################
## Load atlas ##
################

# Load SingleCellExperiment
sce_atlas  <- readRDS(io$atlas.sce)

# Rename genes from ENSEMBLE IDs to gene names
gene_metadata <- fread(io$gene_metadata) %>% .[,c("ens_id","symbol")] %>%
  .[ens_id%in%rownames(sce_atlas) & !is.na(symbol) & symbol!=""]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
sce_atlas <- sce_atlas[rownames(sce_atlas)%in%gene_metadata$ens_id,]
rownames(sce_atlas) <- foo[rownames(sce_atlas)]
stopifnot(sum(is.na(rownames(sce_atlas)))==0)
stopifnot(sum(duplicated(rownames(sce_atlas)))==0)

# Load cell metadata
meta_atlas <- fread(io$atlas.metadata) %>%
  .[stripped==FALSE & doublet==FALSE & stage%in%args$atlas_stages]

# Filter
if (isTRUE(args$test)) meta_atlas <- head(meta_atlas,n=1000)
sce_atlas <- sce_atlas[,meta_atlas$cell] 

################
## Load query ##
################

# Load cell metadata
io$metadata <- paste0(io$basedir,"/results/rna/qc/sample_metadata_after_qc.txt.gz")
meta_query <- fread(io$metadata) %>% 
  # .[pass_rnaQC==TRUE & hybrid_call==FALSE & sample%in%args$query_samples]
  .[pass_QC==TRUE & sample%in%args$query_samples]

# Filter
if (isTRUE(args$test)) meta_query <- head(meta_query,n=1000)

# Load SingleCellExperiment
sce_query <- readRDS(io$sce)[,meta_query$cell]

#############
## Prepare ## 
#############

# Filter out non-expressed genes
sce_query <- sce_query[rowMeans(counts(sce_query))>1e-3,]
sce_atlas <- sce_atlas[rowMeans(counts(sce_atlas))>1e-3,]

# Intersect genes
genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

# Load gene markers to be used as HVGs
# marker_genes.dt <- fread(io$atlas.marker_genes)
# marker_genes.dt <- marker_genes.dt[,head(.SD,n=50),by="celltype"]
# marker_genes <- unique(marker_genes.dt$ens_id)
# marker_genes <- marker_genes[marker_genes%in%genes.intersect]
# print(marker_genes.dt[,.N,by="celltype"])
# stopifnot(all(marker_genes%in%rownames(sce_atlas)))

#########
## Map ##
#########

mapping  <- mapWrap(
  atlas_sce = sce_atlas, atlas_meta = meta_atlas,
  map_sce = sce_query, map_meta = meta_query, 
  genes = NULL, npcs = 50, k = 25
)

head(mapping$mapping)

##########
## Save ##
##########

outfile <- sprintf("%s/mapping_mnn_%s.rds",args$outdir,paste(args$query_samples,collapse="-"))
saveRDS(mapping, outfile)
