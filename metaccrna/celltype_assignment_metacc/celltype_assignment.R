here::i_am("metaccrna/celltype_assignment_metacc/celltype_assignment.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$mofa_model <- paste0(io$basedir, "/results/metaccrna/mofa2/model.rds")
io$umap <- paste0(io$basedir, "/results/metaccrna/mofa2/umap.txt.gz")
io$factors <- paste0(io$basedir, "/results/metaccrna/mofa2/mofa_factors_corrected.txt.gz")
io$outdir <- paste0(io$basedir, "/results/metaccrna/celltype_assignment"); dir.create(io$outdir, showWarnings = F)

opts$k <- 15
opts$input_celltype_column <- "celltype.mapped"
opts$output_celltype_column <- "celltype_predicted"

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE | pass_metQC==TRUE | pass_accQC==TRUE] %>%
  setnames(opts$input_celltype_column,"celltype")

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


######################
## Load coordinates ##
######################

factors.mtx <- fread(io$factors) %>% matrix.please
umap.mtx <- fread(io$umap) %>% matrix.please

###########################################################################
## Create neighbourhood graph for cells that are missing cell type label ##
###########################################################################

cat("Creating neighbourhood graph...\n")

df.observed <- cell_metadata.dt[!is.na(celltype),c("cell","celltype")]
cells.missing.label <- cell_metadata.dt[is.na(celltype),cell]

X.observed <- factors.mtx[df.observed$cell,]
X.missing <- factors.mtx[cells.missing.label,]

# returns a list of 2 (N,K) matrices.:
# - nn.idx: 1-indexed indices
# - nn.dists: distances
knnObj <- nabor::knn(data = X.observed, query = X.missing, k = opts$k)

####################
## kNN prediction ##
####################

cat("Doing kNN prediction...\n")

neighbour.celltypes <- apply(knnObj$nn.idx, 1, function(x) df.observed$celltype[x]) %>% t
predicted.celltype <- apply(neighbour.celltypes, 1, function(x) getmode(x, 1:length(x)))

############################
## Update sample metadata ##
############################

cat("Updating sample metadata...\n")

celltype_predictions.dt <- data.table(
  cell = cells.missing.label,
  tmp = predicted.celltype
)

tmp <- fread(io$metadata)
if (opts$output_celltype_column%in%colnames(tmp)) tmp[[args$output_celltype_column]] <- NULL

cell_metadata_updated.dt <- tmp %>% 
  merge(celltype_predictions.dt, by="cell", all.x=TRUE) %>%
  .[is.na(tmp),tmp:=eval(as.name(opts$input_celltype_column))] %>%
  setnames("tmp",opts$output_celltype_column)

# parse metadata
# cell_metadata_updated.dt[pass_rnaQC==TRUE & pass_rnaQC==FALSE & !is.na(celltype.predicted),doublet_call:=FALSE]

# Save metadata
fwrite(cell_metadata_updated.dt, file.path(io$outdir,"sample_metadata_after_celltype_assignment.txt.gz"), sep="\t", na="NA", quote=F)

#################
## Print stats ##
#################

# print(sprintf("Number of cells that passed met QC but did not have cell type assignment (before imputation): N=%s",tmp[pass_metQC==TRUE & pass_rnaQC==FALSE,.N]))
# print(sprintf("Number of cells that pass met QC and now have cell type assignment (after imputation): N=%s",cell_metadata_updated.dt[pass_metQC==TRUE & pass_rnaQC==FALSE,.N]))

##########
## Plot ##
##########

cat("Plotting...\n")
    
to.plot <- umap.mtx %>%
  as.data.table(keep.rownames = T) %>%
  setnames(c("cell","umap1","umap2")) %>%
  merge(cell_metadata_updated.dt, by="cell") %>%
  setnames(opts$input_celltype_column,"celltype") %>%
  setnames(opts$output_celltype_column,"celltype_predicted")

p1 <- ggplot(to.plot[!is.na(celltype)], aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=celltype), size=1.25, shape=21, color="black", stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  labs(title=sprintf("Cells that have cell type assignment (%s, N=%s)",opts$input_celltype_column,to.plot[!is.na(celltype),.N])) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5, size=rel(0.9)))

p2 <- ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=is.na(celltype), size=is.na(celltype), alpha=is.na(celltype)), shape=21, color="black", stroke=0.05) +
  scale_size_manual(values=c("TRUE"=0.75, "FALSE"=0.4)) +
  scale_fill_manual(values=c("TRUE"="red", "FALSE"="gray60")) +
  scale_alpha_manual(values=c("TRUE"=0.75, "FALSE"=0.25)) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  labs(title=sprintf("Highlighting cells that do not have cell type assignment (in red, N=%s)",to.plot[is.na(celltype),.N])) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5, size=rel(0.9)))
  
p3 <- ggplot(to.plot[is.na(celltype)], aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=celltype_predicted), size=1.25, shape=21, color="black", stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  labs(title=sprintf("Subset of cells after celltype prediction (%s, N=%s)",opts$output_celltype_column,to.plot[is.na(celltype),.N])) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5, size=rel(0.9)))

p <- cowplot::plot_grid(plotlist=list(p1,p2,p3), nrow = 1)

pdf(file.path(argios$outdir,"celltype_assignment.pdf"), width=16, height=6)
print(p)
dev.off()


