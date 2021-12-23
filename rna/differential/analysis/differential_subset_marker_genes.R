#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$differential_results <- file.path(io$basedir,"results_new/rna/differential")
io$outdir <- file.path(io$basedir,"results_new/rna/differential/marker_genes"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

opts$celltypes <- c(
  "Primitive_Streak",
  "Nascent_mesoderm",
  "Mesenchyme",
  "ExE_mesoderm",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors",
  # "Surface_ectoderm",
  "early_Erythroid"
)

# Significance thresholds
opts$min.log2FC <- 1
opts$fdr <- 0.01

##########################################
## Load differential expression results ##
##########################################

# i <- "Gut"; j <- "NMP"
diff.dt <- opts$celltypes %>% map(function(i) {
  file <- file.path(io$differential_results,sprintf("%s_WT_vs_KO.txt.gz",i))
  fread(file) %>% 
    .[abs(logFC)>=opts$min.log2FC & padj_fdr<=opts$fdr] %>%
    .[,celltype:=i] %>%
    return
  }) %>% rbindlist %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$atlas.marker_genes)
genes_to_use <- genes.intersect[genes.intersect%in%unique(marker_genes.dt$gene)]

diff_markers.dt <- diff.dt %>%
  .[gene%in%unique(marker_genes.dt$gene)] %>%
  setorderv("padj_fdr")

# Save
fwrite(diff_markers.dt, file.path(io$outdir,"differential_marker_genes.tsv.gz"), quote=F, sep="\t", na="NA")

#########################################################
## Plot number of diff expr marker genes per cell type ##
#########################################################

to.plot <- diff_markers.dt %>% .[,.N,by="celltype"]

p <- ggbarplot(to.plot, x="celltype", y="N", fill="celltype") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of marker genes") +
  theme(
    axis.text.y = element_text(size=rel(0.65)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    axis.title = element_text(colour="black",size=rel(0.75)),
    axis.ticks.x = element_blank(),
    legend.position = "none"
)

pdf(file.path(io$outdir,"barplot_number_marker_genes.pdf"), width = 8, height = 4)
print(p)
dev.off()

