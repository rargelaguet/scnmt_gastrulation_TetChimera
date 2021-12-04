suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--counts',        type="character",                    help='Counts file')
p$add_argument('--outdir',       type="character",                    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

here::i_am("rna/processing/1_create_seurat_rna.R")
source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
args <- list()
args$counts <- file.path(io$basedir,"processed/rna/counts_GRCm38_100.tsv.gz")
args$metadata <- file.path(io$basedir,"processed/sample_metadata_sce.txt")
args$outdir <- file.path(io$basedir,"processed/rna")
## END TEST ##

#######################
## Load count matrix ##
#######################

rna_counts.mtx <- fread(args$count) %>% matrix.please

##########################
## Load sample metadata ##
##########################

metadata <- fread(args$metadata)

mean(colnames(rna_counts.mtx)%in%metadata$id_rna)
mean(colnames(rna_counts.mtx)%in%metadata$cell)

##################
## Filter genes ##
##################

# Keep protein-coding genes
# if (!is.null(opts$subset.proteincoding)){
#     genes <- fread(opts$subset.proteincoding)[,ens_id]
#     genes <- genes[genes %in% mouse.genes]
#     mouse.genes <- mouse.genes[mouse.genes %in% genes]
#     rna_counts.mtx <- rna_counts.mtx[mouse.genes,]
# }

# Remove duplicated genes
rna_counts.mtx <- rna_counts.mtx[!duplicated(rownames(rna_counts.mtx)),]

# Sanity checks
stopifnot(sum(duplicated(rownames(rna_counts.mtx)))==0)
stopifnot(sum(duplicated(colnames(rna_counts.mtx)))==0)

##########################
## Create Seurat object ##
##########################

cell.info.to.seurat <- cell.info[cell%in%colnames(rna_counts.mtx)] %>% setkey(cell) %>% .[colnames(rna_counts.mtx)] %>% as.data.frame
rownames(cell.info.to.seurat) <- cell.info.to.seurat$cell
stopifnot(rownames(cell.info.to.seurat)==colnames(rna_counts.mtx))
stopifnot(sum(is.na(rownames(cell.info.to.seurat$cell)))==0)

seurat <- CreateSeuratObject(rna_counts.mtx, meta.data = cell.info.to.seurat)

head(seurat@meta.data)

# Add mitochondrial percenatge
seurat[["mitochondrial_percent_RNA"]] <- PercentageFeatureSet(seurat, pattern = "^mt-") %>% round(2)

# Add ribosomal RNA content
ribo.genes <- grep(pattern = "^Rp[l|s]", x = rownames(seurat), value = TRUE)
seurat[["ribosomal_percent_RNA"]] <- PercentageFeatureSet(seurat, features = ribo.genes) %>% round(2)

#####################
## Create metadata ##
#####################

metadata <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL] %>%
  .[,c("cell","barcode","sample","nFeature_RNA","nCount_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA")]

# metadata[,stage:=substr(sample,1,4)]
metadata[,stage:=strsplit(sample,"_") %>% map_chr(1)]

##########
## Save ##
##########

fwrite(metadata, file.path(args$outdir,"metadata.txt.gz"), quote=F, na="NA", sep="\t")
# fwrite(cell.info, paste0(args$outdir,"/cell_info.txt.gz"), quote=F, na="NA", sep="\t")
# fwrite(gene.info, paste0(args$outdir,"/gene_info.txt.gz"), quote=F, na="NA", sep="\t")
saveRDS(seurat, file.path(args$outdir,"seurat.rds"))

