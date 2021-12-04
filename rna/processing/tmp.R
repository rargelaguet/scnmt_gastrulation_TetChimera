suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))

here::i_am("rna/processing/1_create_seurat_rna.R")
source(here::here("settings.R"))
source(here::here("utils.R"))

args <- list()

## END TEST ##

###############
## Load data ##
###############

sce <- readRDS(file.path(io$basedir,"processed/rna/SingleCellExperiment.rds"))

# metadata <- colData(sce) %>% as.data.table(keep.rownames = T) %>% setnames("rn","id_rna")

gene_metadata.dt <- fread(io$gene_metadata) %>% 
  .[!is.na(symbol) & symbol!=""] %>% 
  .[!duplicated(symbol)]

##################
## Rename genes ##
##################

# Filter common genes
genes <- intersect(rownames(sce),gene_metadata.dt$ens_id)
sce <- sce[genes,]
gene_metadata.dt <- gene_metadata.dt[ens_id%in%genes] %>% setkey(ens_id) %>% .[genes]

stopifnot(rownames(sce)==gene_metadata.dt$ens_id)

tmp <- gene_metadata.dt$symbol; names(tmp) <- gene_metadata.dt$ens_id
rownames(sce) <- tmp[rownames(sce)]

stopifnot(!is.na(rownames(sce)))
stopifnot(sum(duplicated(rownames(sce)))==0)

##########################
## Create Seurat object ##
##########################

seurat <- as.Seurat(sce)
seurat <- RenameAssays(seurat, "originalexp"="RNA")
head(seurat@meta.data)

# Add RNA stats
# seurat[["nFeature_RNA"]] <- colSums(GetAssayData(seurat,'counts'))
# seurat[["nCount_RNA"]] <- colSums(GetAssayData(seurat,'counts'))

# Add mitochondrial percenatge
seurat[["mit_percent_RNA"]] <- PercentageFeatureSet(seurat, pattern = "^mt-") %>% round(2)

# Add ribosomal RNA content
ribo.genes <- grep(pattern = "^Rp[l|s]", x = rownames(seurat), value = TRUE)
seurat[["rib_percent_RNA"]] <- PercentageFeatureSet(seurat, features = ribo.genes) %>% round(2)

##########
## Save ##
##########

metadata <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL] %>%
  .[,c("nCount_originalexp","nFeature_originalexp","sierra_rna","sierra_dna","i7","i5"):=NULL]

args$outdir <- file.path(io$basedir,"processed/rna_new")
fwrite(metadata, file.path(args$outdir,"metadata.txt.gz"), quote=F, na="NA", sep="\t")
saveRDS(seurat, file.path(args$outdir,"seurat.rds"))
