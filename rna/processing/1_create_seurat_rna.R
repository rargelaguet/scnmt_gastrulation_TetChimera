suppressPackageStartupMessages(library(Seurat))

here::i_am("rna/processing/1_create_seurat_rna.R")
source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                    help='Cell metadata file')
p$add_argument('--counts',        type="character",                    help='Counts file')
p$add_argument('--gene_metadata',        type="character",                    help='Gene metadata file')
p$add_argument('--outdir',       type="character",                    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$counts <- file.path(io$basedir,"processed/rna/counts.tsv.gz")
# args$metadata <- file.path(io$basedir,"processed/sample_metadata.txt.gz")
# args$gene_metadata <- io$gene_metadata
# args$outdir <- file.path(io$basedir,"processed/rna")
## END TEST ##

#######################
## Load count matrix ##
#######################

rna_counts.mtx <- fread(args$count) %>% matrix.please

##########################
## Load sample metadata ##
##########################

metadata <- fread(args$metadata, select=c(1,3,4,5,6,7,8,9,10,11,12,13,14))

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(args$gene_metadata)[,c("ens_id","symbol")] %>% 
  .[ens_id%in%rownames(rna_counts.mtx) & symbol!=""] %>% .[!duplicated(symbol)]

########################################
## Rename genes from symbol to ens_id ##
########################################

genes <- intersect(gene_metadata.dt$ens_id,rownames(rna_counts.mtx))
gene_metadata.dt <- gene_metadata.dt %>% .[ens_id%in%genes] %>% setkey(ens_id) %>% .[genes]
rna_counts.mtx <- rna_counts.mtx[genes,] 
tmp <- gene_metadata.dt$symbol; names(tmp) <- gene_metadata.dt$ens_id
rownames(rna_counts.mtx) <- tmp[rownames(rna_counts.mtx)]

# Sanity checks
stopifnot(!is.na(rownames(rna_counts.mtx)))
stopifnot(!duplicated(rownames(rna_counts.mtx)))

###########################
## Parse sample metadata ##
###########################

# Parse column names
# colnames(metadata) <- make.names(colnames(metadata))

# create sample column
plate2sample <- c(
  "E7.5_tet_chimera_plate3" = "E7.5_TET_TKO",
  "E7.5_tet_crispr_plate5" = "E7.5_TET_TKO_crispr",
  "E7.5_tet_crispr_plate6" = "E7.5_TET_TKO_crispr",
  "E8.5_oct20_plate1" = "E8.5_WT_CD41+",
  "E8.5_oct20_plate2" = "E8.5_TET_TKO_CD41+",
  "E8.5_oct20_plate3" = "E8.5_WT_KDR+",
  "E8.5_oct20_plate4" = "E8.5_TET_TKO_KDR+",
  "E8.5_oct20_plate5" = "E8.5_WT_KDR+_CD41+",
  "E8.5_oct20_plate6" = "E8.5_TET_TKO_KDR+_CD41+",
  "E8.5_oct20_plate7" = "E8.5_WT",
  "E8.5_oct20_plate8" = "E8.5_TET_TKO"
)
metadata[,sample:=stringr::str_replace_all(plate,plate2sample)]
table(metadata$sample)

# Remove id_acc for samples that were processed with MT-seq
metadata[method=="mt",id_acc:=NA]

# Remove unused columns
metadata[,c("tdTOM","KDR-Cy7","CD41-BV421"):=NULL]

# Sanity checks
stopifnot(sort(colnames(rna_counts.mtx))==sort(metadata$id_rna))
# all(metadata[!is.na(id_acc),id_acc]%in%gsub(".tsv.gz","",list.files(io$acc_data_raw, pattern = "*.tsv.gz")))
# all(metadata[!is.na(id_met),id_met]%in%gsub(".tsv.gz","",list.files(io$met_data_raw, pattern = "*.tsv.gz")))

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

metadata_to_seurat <- metadata %>% setkey(id_rna) %>% .[colnames(rna_counts.mtx)] %>% as.data.frame
rownames(metadata_to_seurat) <- metadata_to_seurat$id_rna
stopifnot(rownames(metadata_to_seurat)==colnames(rna_counts.mtx))

seurat <- CreateSeuratObject(rna_counts.mtx, meta.data = metadata_to_seurat)

# head(seurat@meta.data)

# Add mitochondrial percenatge
seurat[["mit_percent_RNA"]] <- PercentageFeatureSet(seurat, pattern = "^mt-") %>% round(2)

# Add ribosomal RNA content
ribo.genes <- grep(pattern = "^Rp[l|s]", x = rownames(seurat), value = TRUE)
seurat[["rib_percent_RNA"]] <- PercentageFeatureSet(seurat, features = ribo.genes) %>% round(2)

#####################
## Create metadata ##
#####################

metadata.to.save <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL]
  
##########
## Save ##
##########

fwrite(metadata.to.save, file.path(args$outdir,"metadata.txt.gz"), quote=F, na="NA", sep="\t")
saveRDS(seurat, file.path(args$outdir,"seurat.rds"))

