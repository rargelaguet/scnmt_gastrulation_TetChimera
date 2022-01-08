source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results_new/rna/individual_genes/pseudobulk")
io$sce.pseudobulk <- file.path(io$basedir,"results_new/rna/pseudobulk/SingleCellExperiment_pseudobulk_ko_celltype.rds")

# Define options
opts$celltypes = c(
  # "Spinal_cord",
  # "Surface_ectoderm",
  # "Gut",
  # "Pharyngeal_mesoderm",
  # "Endothelium",
  "Haematoendothelial_progenitors",
  "Blood_progenitors",
  "early_Erythroid",
  "late_Erythroid"
)
opts$min.cells <- 10

# opts$rename <- c(
#   "Haematoendothelial_progenitors" = "Haemato._progenitors"
# )

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$class <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)

# Filter
sce <- sce[,sce$celltype%in%opts$celltypes]

# Filter by minimum number of cells
# sce <- sce[,names(which(metadata(sce)$n_cells>=opts$min.cells))]

metadata.dt <- data.table(
  sample = colnames(sce),
  celltype = factor(sce$celltype, levels=opts$celltypes),
  class = sce$class
)

###############
## Lineplots ##
###############

genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")

to.plot <- logcounts(sce)[genes.to.plot,]  %>%
  reshape2::melt() %>% as.data.table() %>%
  setnames(c("gene","sample","expr")) %>%
  .[,gene_class:=ifelse(grepl("Dnmt",gene),"DNMTs","TETs")] %>%
  merge(metadata.dt,by="sample") 

to.plot %>%
  # .[,celltype:=stringr::str_replace_all(celltype,opts$rename)] %>%
  .[,tmp:=sprintf("%s (%s)",gene_class,class)] 

ggline(to.plot, x="celltype", y="expr", color="gene") +
  # geom_bar(stat="identity", color="black") +
  scale_color_brewer(palette="Dark2") +
  facet_wrap(~tmp, scales="fixed") +
  labs(x="",y="Gene expression") +
  theme_classic() +
  guides(x = guide_axis(angle = 90)) +
  ggrepel::geom_text_repel(aes_string(label="gene"), nudge_x=0.25, data=to.plot[celltype=="late_Erythroid"]) +
  # geom_text(aes_string(label="gene"), position=position_nudge(x=0.25), data=to.plot[celltype=="late_Erythroid"]) +
  theme(
    strip.text = element_text(size=rel(0.85)),
    axis.text.x = element_text(colour="black",size=rel(0.9)),
    axis.text.y = element_text(colour="black",size=rel(0.9)),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(colour="black",size=rel(1.0)),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size=rel(0.85))
  )

# pdf(outfile, width=10, height=9)
# print(p)
# dev.off()

