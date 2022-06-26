# here::i_am("rna_atac/mofa/plot_mofa_results.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$mofa_model <- paste0(io$basedir, "/results/metaccrna/mofa2/model.rds")
io$outdir <- paste0(io$basedir, "/results/metaccrna/mofa2"); dir.create(io$outdir, showWarnings = F)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE | pass_metQC==TRUE | pass_accQC==TRUE] %>%
  setnames("celltype.mapped","celltype")

###############
## Load MOFA ##
###############

if (tools::file_ext(io$mofa_model)=="rds") {
  MOFAobject <- readRDS(io$mofa_model)
} else if (file_ext(io$mofa_model)=="hdf5") {
  MOFAobject <- load_model(io$mofa_model, load_data = F)
}

# Add sample metadata
cells <- as.character(unname(unlist(MOFA2::samples_names(MOFAobject))))
sample_metadata_to_mofa <- copy(cell_metadata.dt) %>%
  setnames("cell","sample") %>% .[sample%in%cells] %>% setkey(sample) %>% .[cells]
stopifnot(all(cells==sample_metadata_to_mofa$cell))
samples_metadata(MOFAobject) <- sample_metadata_to_mofa

table(MOFAobject@samples_metadata$celltype)

#############################
## Plot variance explained ##
#############################

p <- plot_variance_explained(MOFAobject, plot_total = T)[[2]]

pdf(sprintf("%s/mofa_var_explained_total.pdf",io$outdir), width=6, height=3)
print(p)
dev.off()

p <- plot_variance_explained(MOFAobject, factors = 1:MOFAobject@dimensions$K, x="view", y="factor", max_r2 = 5) +
  theme(legend.position = "top")

pdf(sprintf("%s/mofa_var_explained.pdf",io$outdir), width=6, height=8)
print(p)
dev.off()

##################
## Plot factors ##
##################

p <- plot_factors(MOFAobject, factors = c(8,9), color_by = "celltype") +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none"
  )

pdf(sprintf("%s/mofa_scatterplot_factors_12.pdf",io$outdir), width=6, height=5)
print(p)
dev.off()

##########################################
## Extract factors and batch correction ##
##########################################

library(batchelor)

factors.to.use <- 1:get_dimensions(MOFAobject)[["K"]]
factors.mtx <- get_factors(MOFAobject, factors=factors.to.use)[[1]]

opts$batch_correction <- TRUE

if (opts$batch_correction) {
  batch_variable <- MOFAobject@samples_metadata$method
  factors_corrected.mtx <- reducedMNN(factors.mtx, batch=batch_variable)$corrected 
  # colnames(Z_corrected) <- colnames(Z)
} else {
  factors_corrected.mtx <- factors.mtx
}

factors_corrected.dt <- factors_corrected.mtx %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
fwrite(factors_corrected.dt, file.path(io$outdir,"mofa_factors_corrected.txt.gz"))

##########
## UMAP ##
##########

# Run
umap.mtx <- uwot::umap(factors_corrected.mtx, n_neighbors=25, min_dist=0.50, metric="cosine")
umap.dt <- umap.mtx %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
fwrite(umap.dt, file.path(io$outdir,"umap.txt.gz"))

# Plot
to.plot <- umap.dt %>%
  .[,sample:=rownames(factors.mtx)] %>%
  merge(MOFAobject@samples_metadata %>% as.data.table)

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=celltype)) +
  # geom_point(size=1, shape=21, stroke=0.1) +
  ggrastr::geom_point_rast(size=1, shape=21, stroke=0.1) +
  scale_fill_manual(values=opts$celltype.colors) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.position="none"
  )

pdf(file.path(io$outdir,"mofa_umap_celltype.pdf"), width=7, height=5)
print(p)
dev.off()


##########
## TEST ##
##########

# ggplot(to.plot, aes(x=V1, y=V2, fill=method)) +
#   geom_point(size=1, shape=21, stroke=0.1) +
#   # scale_fill_manual(values=opts$celltype.colors) +
#   guides(fill = guide_legend(override.aes = list(size=2))) +
#   theme_classic() +
#   ggplot_theme_NoAxes() +
#   theme(
#     legend.position="right"
#   )
