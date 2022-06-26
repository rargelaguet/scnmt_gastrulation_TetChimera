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