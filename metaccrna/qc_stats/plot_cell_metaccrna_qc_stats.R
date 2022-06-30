library(VennDiagram)

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/metaccrna/qc_stats"); dir.create(io$outdir, showWarnings = F)

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata) %>% .[ko_type!="crispr"]

##################
## Venn Diagram ##
##################

to.plot <- cell_metadata.dt %>%
  .[pass_rnaQC==T | pass_metQC==T | pass_accQC==T]

p <- venn.diagram(
  x = list("RNA expression"=to.plot[pass_rnaQC==T,cell],
           "DNA methylation"=to.plot[pass_metQC==T,cell],
           "Chr. accessibility"=to.plot[pass_accQC==T,cell]
           ),
  filename=NULL,
  col="transparent", fill=c("#3CB54E","#F8766D","#00BFC4"), alpha = 0.60, cex = 1.5, 
  fontfamily = "serif", fontface = "bold")

pdf(file.path(io$outdir,"venn_metaccrna_cell_numbers.pdf"))
grid.draw(p)
dev.off()

#############################################
## Correlate RNA and epigenetic QC metrics ##
#############################################

to.plot <- cell_metadata.dt#[,c("cell","method","nFeature_RNA","nCount_RNA","nCG","nGC")]

p <- ggscatter(to.plot, x="nCG", y="nFeature_RNA", color="pass_rnaQC",
               add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
  # scale_fill_manual(values=opts$celltype.colors) +
  stat_cor(method = "pearson") +
  # facet_wrap(~class, scales="fixed") +
  theme_classic() +
  labs(y="", x="") +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=rel(0.90))
  )

pdf(file.path(io$outdir,"scatterplot_diff_proportions.pdf"), width=8.5, height=5)
print(p)
dev.off()
