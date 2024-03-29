here::i_am("plot_individual_genes/pseudobulk/plot_individual_genes_pseudobulk.R")

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
  "Spinal_cord",
  "Surface_ectoderm",
  "Gut",
  "Pharyngeal_mesoderm",
  "Endothelium",
  "Haematoendothelial_progenitors",
  "Blood_progenitors",
  "early_Erythroid",
  "late_Erythroid"
)
opts$min.cells <- 10

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$class <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)

# Filter celltypes
# sce <- sce[,sce$celltype%in%opts$celltypes]

# Filter by minimum number of cells
sce <- sce[,names(which(metadata(sce)$n_cells>=opts$min.cells))]

#############
## Barplots ##
#############

# genes.to.plot <- c("Tex19.1","Morc1","Dppa3","Rex1","Dppa5a","Dppa4","Dppa2","Zfp981")
# genes.to.plot <- c("Pou5f1","Epcam","Fst","Lefty2","Cdkn1c","Acta1")
# genes.to.plot <- grep("^Hox",rownames(sce),value=T)
genes.to.plot <- fread(io$atlas.marker_genes)[,gene] %>% unique

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  if (gene %in% rownames(sce)) {
    print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    outfile <- sprintf("%s/%s_barplot_pseudobulk.pdf",io$outdir,gene)
  
    to.plot <- data.table(
      sample = colnames(sce),
      expr = logcounts(sce)[gene,],
      class = sce$class,
      celltype = sce$celltype
    )
    
    p <- ggplot(to.plot, aes(x=class, y=expr, fill=class)) +
      geom_bar(stat="identity", color="black") +
      scale_fill_brewer(palette="Dark2") +
      facet_wrap(~celltype, scales="fixed") +
      theme_classic() +
      labs(x="",y=sprintf("%s expression",gene)) +
      theme(
        strip.text = element_text(size=rel(0.85)),
        axis.text.x = element_text(colour="black",size=rel(0.9)),
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=rel(1.0)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.85))
      )
    
      pdf(outfile, width=10, height=9)
      print(p)
      dev.off()
        
  } else {
    print(sprintf("%s not found",gene))
  }
}

