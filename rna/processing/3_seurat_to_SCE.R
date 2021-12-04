suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--test',            action="store_true",                 help='Testing mode')
p$add_argument('--normalise',       action="store_true",                 help='Log-Normalise?')
p$add_argument('--seurat',         type="character", help='Seurat object (input)')
p$add_argument('--metadata',         type="character", help='Metadata file')
p$add_argument('--outfile',         type="character", help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))

## START TEST ##
# args <- list()
# args$outfile <- file.path(io$basedir,"processed/rna_new/SingleCellExperiment.rds")
# args$metadata <- file.path(io$basedir,"results/rna_new/qc/sample_metadata_after_qc.txt.gz")
# args$seurat <- file.path(io$basedir,"processed/rna_new/seurat.rds")
# args$normalise <- FALSE
# args$test <- FALSE
## END TEST ##

# Sanity checks
stopifnot(args$samples%in%opts$samples)
if (args$test) args$samples <- head(args$samples,n=2)

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(args$metadata) %>% .[pass_rnaQC==TRUE]
table(sample_metadata$sample)

# Load seurat
seurat <- readRDS(args$seurat)[,sample_metadata$cell]

#####################################
## Convert to SingleCellExperiment ##
#####################################

sce <- as.SingleCellExperiment(seurat)

# remove logcounts assays
sce@assays@data[["logcounts"]] <- NULL

# Add metadata
# stopifnot(sample_metadata$cell%in%colnames(sce))
# stopifnot(colnames(sce)%in%sample_metadata$cell)
sample_metadata <- sample_metadata %>% .[cell%in%colnames(sce)] %>% setkey(cell) %>% .[colnames(sce)]
stopifnot(sample_metadata$cell == colnames(sce))
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

##########################
## Compute size factors ##
##########################

clusts <- as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))
# clusts <- as.numeric(quickCluster(sce))
min.clust <- min(table(clusts))/2
new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce <- computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

###################
## Log Normalise ##
###################

if (args$normalise) {
	sce <- logNormCounts(sce)
}

##########
## Plot ##
##########

# to.plot <- data.frame(X = Matrix::colSums(counts(sce)), Y = sizeFactors(sce))
# ggplot(to.plot, mapping = aes(x = X, y = Y)) +
#   geom_point() +
#   labs(x = "Number of UMIs", y = "Size Factor") +
#   theme_classic()

##########
## Save ##
##########

saveRDS(sce, args$outfile)