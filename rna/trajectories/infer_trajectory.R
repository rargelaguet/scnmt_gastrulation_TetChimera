here::i_am("rna/trajectories/infer_trajectory.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(destiny))
suppressPackageStartupMessages(library(reticulate))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='') 
p$add_argument('--sce',  type="character",              help='') 
p$add_argument('--trajectory_name',  type="character",              help='') 
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
args <- list()
args$sce <- io$rna.sce
args$metadata <- io$metadata
args$trajectory_name <- "blood"
args$outdir <- file.path(io$basedir,"results_new/rna/trajectories/blood")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

# Options
opts$celltype_trajectory_dic <- list(
  "blood" = c("Haematoendothelial_progenitors", "Blood_progenitors", "early_Erythroid", "late_Erythroid"),
  "ectoderm" = c("Epiblast", "Rostral_neurectoderm", "Forebrain_Midbrain_Hindbrain"),
  "endoderm" = c("Epiblast", "Anterior_Primitive_Streak", "Def._endoderm", "Gut"),
  "mesoderm" = c("Epiblast", "Primitive_Streak", "Nascent_mesoderm")
)

stopifnot(args$trajectory_name%in%names(opts$celltype_trajectory_dic))
opts$celltypes <- opts$celltype_trajectory_dic[[args$trajectory_name]]

opts$genes2plot <- list(
  "blood" = c("Hbb-y","Etv2"),
  "ectoderm" = c("Utf1","Crabp1"),
  "endoderm" = c("Pou5f1","Krt8"),
  "mesoderm" = c("Dnmt3b","Mesp1")
)
stopifnot(args$trajectory_name%in%names(opts$genes2plot))
opts$genes_to_plot <- opts$genes2plot[[args$trajectory_name]]

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  # .[grep("WT",class)] %>%
  setnames("celltype.mapped","celltype") %>%
  .[pass_rnaQC==TRUE & celltype%in%opts$celltypes] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

table(sample_metadata$stage)
table(sample_metadata$celltype)

#########################
## Load RNA expression ##
#########################
  
# Load SingleCellExperiment
sce <- load_SingleCellExperiment(
  file = args$sce, 
  normalise = TRUE, 
  cells = sample_metadata$id_rna, 
  remove_non_expressed_genes = TRUE
)

colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("id_rna") %>%
  .[colnames(sce),] %>% DataFrame()

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- rownames(decomp)[decomp$p.value <= 0.01]

# Subset HVGs
sce_filt <- sce[hvgs,]
dim(sce_filt)

#########
## PCA ##
#########

sce_filt <- runPCA(sce_filt, ncomponents = 5, ntop=nrow(sce_filt))

# Plot PCA
pdf(file.path(args$outdir,sprintf("%s_pca_celltype.pdf",args$trajectory_name)), width=8, height=5)
plotPCA(sce_filt, colour_by="celltype", ncomponents = c(1,2)) +
  scale_colour_manual(values=opts$celltype.colors) + theme(legend.position="none")
dev.off()


###################
## Diffusion map ##
###################

set.seed(42)
dm <- DiffusionMap(sce_filt, n_pcs=2)

# Add to the SingleCellExperiment object
reducedDim(sce_filt, "DiffusionMap") <- dm@eigenvectors[,c(1,2)]

# Plot
pdf(file.path(args$outdir,sprintf("%s_diffmap_celltype.pdf",args$trajectory_name)), width=8, height=5)
plotReducedDim(sce_filt, dimred = "DiffusionMap", colour_by="celltype", ncomponents = c(1,2)) +
  scale_colour_manual(values=opts$celltype.colors) + theme(legend.position="none")
dev.off()


##########
## UMAP ##
##########

# set.seed(42)
# dm <- RunUMAP(sce_filt, n_pcs=2)

# Add to the SingleCellExperiment object
# reducedDim(sce_filt, "DiffusionMap") <- dm@eigenvectors[,c(1,2)]

# Plot
# pdf(file.path(args$outdir,sprintf("%s_diffmap_celltype.pdf",args$trajectory_name)), width=8, height=5)
# plotReducedDim(sce_filt, dimred = "DiffusionMap", colour_by="celltype", ncomponents = c(1,2)) +
#   scale_colour_manual(values=opts$celltype.colors) + theme(legend.position="none")
# dev.off()

####################
## Prepare output ##
####################

pseudotime.dt <- data.table(
  id_rna = colnames(sce_filt),
  PC1 = reducedDim(sce_filt,"PCA")[,"PC1"],
  PC2 = reducedDim(sce_filt,"PCA")[,"PC2"],
  DC1 = reducedDim(sce_filt,"DiffusionMap")[,"DC1"],
  DC2 = reducedDim(sce_filt,"DiffusionMap")[,"DC2"]
)

#########################################################################
## Scatterplots of pseudotime values versus expression of marker genes ##
#########################################################################

# Denoise
# pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[colnames(sce),]
# assay(sce_filt,"logcounts_denoised") <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce_filt)), D=pdist(pca.rna), k=25)

rna.dt <- data.table(as.matrix(logcounts(sce)[opts$genes_to_plot,]), keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="id_rna", value.name="expr") %>%
  merge(sample_metadata[,c("id_rna","celltype")])

to.plot <- pseudotime.dt %>%
  merge(rna.dt[gene%in%opts$genes_to_plot],by="id_rna") %>%
  .[,pca:=minmax.normalisation(PC1)] %>% .[,diffmap:=minmax.normalisation(DC1)] %>%
  .[,gene:=factor(gene,levels=opts$genes_to_plot)] %>%
  .[,expr:=minmax.normalisation(expr),by="gene"]

for (i in c("pca","diffmap")) {
  
  p <- ggplot(to.plot, aes_string(x=i, y="expr")) +
    # geom_point(aes(fill=celltype), size=1.75, shape=21, stroke=0.1, alpha=0.5) +
    ggrastr::geom_point_rast(aes(fill=celltype), size=2.5, shape=21, stroke=0.1, alpha=0.5) +
    stat_smooth(aes(fill=expr), method="loess", color="black", alpha=0.85, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    # viridis::scale_fill_viridis() +
    scale_color_manual(values=opts$celltype.colors) +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~gene, scales="free_y") +
    labs(x=sprintf("Pseudotime (%s)",i), y="Gene expression") +
    theme_classic() +
    guides(fill="none", color="none") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      # axis.text.y = element_text(color="black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_text(size=rel(1.25), color="black"),
      legend.position="right",
      legend.title = element_blank()
    )
  
  pdf(file.path(args$outdir,sprintf("%s_%s_vs_expression.pdf",args$trajectory_name,i)), width=12, height=8)
  print(p)
  dev.off()
}

#################
## Save output ##
#################

# Save pseudotime
fwrite(pseudotime.dt, file.path(args$outdir,sprintf("/%s_trajectory.txt.gz",args$trajectory_name)))

# Save metadata
# tmp <- sample_metadata %>% merge(pseudotime.dt,by="id_rna",all.x=TRUE)
fwrite(sample_metadata, file.path(args$outdir,sprintf("%s_sample_metadata.txt.gz",args$trajectory_name)), sep="\t", na="NA")

# Save SingleCellExperiment
trajectory.sce <- sce
reducedDims(trajectory.sce) <- reducedDims(sce_filt)
logcounts(trajectory.sce) <- NULL
saveRDS(trajectory.sce, file.path(args$outdir,sprintf("%s_SingleCellExperiment.rds",args$trajectory_name)))

