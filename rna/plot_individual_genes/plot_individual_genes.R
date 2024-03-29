here::i_am("rna/plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O ##
io$outdir <- file.path(io$basedir,"results/rna/individual_genes"); dir.create(io$outdir, showWarnings = F)

## Define options ##

# Define cell types to plot
# opts$celltypes = c(
#   "Epiblast",
#   "Primitive_Streak",
#   "Caudal_epiblast",
#   "PGC",
#   "Anterior_Primitive_Streak",
#   "Notochord",
#   "Def._endoderm",
#   "Gut",
#   "Nascent_mesoderm",
#   "Mixed_mesoderm",
#   "Intermediate_mesoderm",
#   "Caudal_Mesoderm",
#   "Paraxial_mesoderm",
#   "Somitic_mesoderm",
#   "Pharyngeal_mesoderm",
#   "Cardiomyocytes",
#   "Allantois",
#   "ExE_mesoderm",
#   "Mesenchyme",
#   "Haematoendothelial_progenitors",
#   "Endothelium",
#   "Blood_progenitors",
#   "early_Erythroid",
#   "late_Erythroid",
#   "NMP",
#   "Neurectoderm",
#   "Neural_crest",
#   "Forebrain_Midbrain_Hindbrain",
#   "Spinal_cord",
#   "Surface_ectoderm"
# )

opts$celltypes = c(
  # "Pharyngeal_mesoderm",
  # "ExE_mesoderm",
  # "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors",
  "early_Erythroid"
  # "late_Erythroid",
  # "Surface_ectoderm"
)

# Define samples to plot
opts$samples <- c(
  "E7.5_TET_TKO",
  # "E7.5_TET_TKO_crispr",
  "E8.5_WT_CD41+",
  "E8.5_TET_TKO_CD41+",
  "E8.5_WT_KDR+",
  "E8.5_TET_TKO_KDR+",
  "E8.5_WT_KDR+_CD41+",
  "E8.5_TET_TKO_KDR+_CD41+",
  "E8.5_WT",
  "E8.5_TET_TKO"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes & sample%in%opts$samples] %>%
  .[,sample:=factor(sample,levels=opts$samples)] %>%
  .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes)]

sample_metadata %>% 
  .[,ko:=ifelse(grepl("KO",sample),"Tet-TKO","WT")] %>%
  .[ko_type=="crispr",ko:="Tet-TKO (crispr)"]

table(sample_metadata$ko)
table(sample_metadata$celltype.mapped)

# Only consider cell types with sufficient observations in both KO and WT cells
celltypes.to.plot <- sample_metadata[,.(N=.N),by=c("ko","celltype.mapped")] %>% .[N>=10] %>% .[,.N,by="celltype.mapped"] %>% .[N>1,celltype.mapped]
sample_metadata <- sample_metadata[celltype.mapped%in%celltypes.to.plot]

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$rna.sce, cells=sample_metadata$id_rna, normalise = TRUE)

# Add sample metadata as colData
# colData(sce) <- sample_metadata %>% tibble::column_to_rownames("id_rna") %>% DataFrame

##########
## Plot ##
##########

genes.to.plot <- c("Cited4","Runx1","Klf1")
# genes.to.plot <- fread("/Users/argelagr/data/tet_chimera_nmtseq/results_new/rna/differential/marker_genes/differential_marker_genes.tsv.gz")$gene %>% unique
# genes.to.plot <- c("Lefty1","Cd34","Tmsb4x","Fgf3","Spata7","Cer1","Spink1","Dppa4","Dppa5a","Prc1","Lefty2","Ube2c","Hba-x","Hbb-y","Hba-a1","Hbb-bh1")
# genes.to.plot <- c("Vegfa","Vegfb","Vegfc","Vegfd","Kdr","Flt1","Tal1","Runx1","Etv2)
# genes.to.plot <- c("Tet1","Tet2","Tet3","Dnmt1","Dnmt3a","Dnmt3b","Dnmt3l")
# genes.to.plot <- rownames(sce)[grep("tomato",rownames(sce))]
# genes.to.plot <- fread(io$atlas.marker_genes)$gene %>% unique %>% .[!grepl("Rik$",.)]
# genes.to.plot <- fread("/Users/ricard/data/gastrulation10x/results/differential/celltypes/E8.5/Neural_crest_vs_Forebrain_Midbrain_Hindbrain.txt.gz") %>% .[sig==T & logFC<(-2),gene]

stopifnot(genes.to.plot%in%rownames(sce))

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  if (gene %in% rownames(sce)) {
    print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    
    to.plot <- data.table(
      id_rna = colnames(sce),
      expr = logcounts(sce)[gene,]
    ) %>% merge(sample_metadata[,c("id_rna","ko","celltype.mapped")], by="id_rna") %>%
      .[,N:=.N,by=c("ko","celltype.mapped")] %>% .[N>=10]
    
    # Plot KO vs WT, facet by cell type
    
    p <- ggplot(to.plot, aes(x=ko, y=expr, fill=ko)) +
      geom_violin(scale = "width", alpha=0.8) +
      geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
      stat_summary(fun.data = give.n, geom = "text", size=3) +
      # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
      scale_fill_manual(values=opts$class.colors) +
      # scale_fill_brewer(palette="Dark2") +
      facet_wrap(~celltype.mapped, scales="fixed") +
      theme_classic() +
      labs(x="",y=sprintf("%s expression",gene)) +
      # guides(x = guide_axis(angle = 90)) +
      theme(
        strip.text = element_text(size=rel(0.85)),
        axis.text.x = element_text(colour="black",size=rel(1)),
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=rel(1.0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(1.0))
      )
    
      pdf(sprintf("%s/%s_per_class.pdf",io$outdir,gene), width=6, height=8)
      print(p)
      dev.off()
      
      # Plot per cell type, facet by WT vs KO 
      
      p <- ggplot(to.plot, aes(x=celltype.mapped, y=expr, fill=celltype.mapped)) +
        geom_violin(scale = "width", alpha=0.8) +
        geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
        stat_summary(fun.data = give.n, geom = "text", size=3) +
        scale_fill_manual(values=opts$celltype.colors) +
        facet_wrap(~ko, nrow = 2, scales="fixed") +
        theme_classic() +
        labs(x="",y=sprintf("%s expression",gene)) +
        guides(x = guide_axis(angle = 90)) +
        theme(
          strip.text = element_text(size=rel(0.85)),
          axis.text.x = element_text(colour="black",size=rel(0.75)),
          axis.text.y = element_text(colour="black",size=rel(0.9)),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(colour="black",size=rel(1.0)),
          legend.position = "none"
        )
      
      pdf(sprintf("%s/%s_per_celltype.pdf",io$outdir,gene), width=8, height=8)
      print(p)
      dev.off()
      
  } else {
    print(sprintf("%s not found",gene))
  }
}

