here::here("metaccrna/plot_individual_features/plot_individual_features.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--rna_gene',           type="character",   help='Feature name for RNA expression (in ENSEMBL ID)')
p$add_argument('--met_feature',    type="character",   help='Feature name for DNA methylation')
p$add_argument('--met_context',    type="character",   help='Genomic context for DNA methylation')
p$add_argument('--acc_feature',    type="character",   help='Feature name for chromatin accessibility')
p$add_argument('--acc_context',    type="character",   help='Genomic context for chromatin accessibility')
p$add_argument('--celltypes',      type="character",   nargs='+',   help='celltypes to plot')
p$add_argument('--outdir',         type="character",              help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$rna_gene <- "Inpp5d"
# args$met_feature <- "chr1:87640323-87640923"
# args$met_context <- "multiome_peaks"
# args$acc_feature <- "chr1:87640323-87640923"
# args$acc_context <- "multiome_peaks"
# args$celltypes <- c("Endothelium","Blood_progenitors","Erythroids")
# args$outdir <- file.path(io$basedir,"results/metaccrna/individual_features")
## END TEST ##

print(args)

#####################
## Define settings ##
#####################

# Define stages to plot

# Define colors for the omics
opts$color <- c(
  "RNA expression"="#3CB54E",
  "Chromatin accessibility"="#00BFC4",
  "DNA methylation"="#F37A71"
)

# Define minimum coverage per cell
opts$min_cpg <- 1
opts$min_gpc <- 3

# Rename cell types
opts$rename.celltypes <- c(
  "early_Erythroid" = "Erythroids",
  "late_Erythroid" = "Erythroids",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak",
  "Mixed_mesoderm" = "Nascent_mesoderm",
  "Allantois" = "ExE_mesoderm"
)

######################
## Load sample data ##
######################

cell_metadata.dt <- fread(io$metadata) %>%
  setnames("celltype.mapped","celltype") %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename.celltypes)] %>%
  .[celltype%in%args$celltypes & !grepl("crispr",class)] %>% .[,celltype:=factor(celltype, levels=args$celltypes)] %>%
  .[,class:=ifelse(grepl("WT",class),"WT","Tet-TKO")] %>%
  .[,celltype_class:=sprintf("%s\n(%s)", celltype,class)]

tmp <- expand.grid(c("WT","Tet-TKO"),args$celltypes)
cell_metadata.dt %>% .[,celltype_class:=factor(celltype_class, levels=sprintf('%s\n(%s)',tmp[,2],tmp[,1]))]

table(cell_metadata.dt$celltype_class)

# Define cells to use
opts$met_cells <- cell_metadata.dt %>% .[pass_metQC==T,id_met]
opts$rna_cells <- cell_metadata.dt %>% .[pass_rnaQC==T,id_rna]
opts$acc_cells <- cell_metadata.dt %>% .[pass_accQC==T,id_acc]

cell_metadata.dt[,c("pass_rnaQC","pass_metQC","pass_accQC"):=NULL]

###############
## Load data ##
###############

# Load DNA methylation data
met.dt <- fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,args$met_context)) %>%
  setnames(c("id_met","id","anno","Nmet","N","value")) %>%
  .[id%in%args$met_feature & id_met%in%opts$met_cells] %>% .[N>=opts$min_cpg]

# Load DNA accessibility data
acc.dt <- fread(sprintf("%s/%s.tsv.gz",io$acc_data_parsed,args$acc_context)) %>%
  setnames(c("id_acc","id","anno","Nmet","N","value")) %>%
  .[id%in%args$acc_feature & id_acc%in%opts$acc_cells] %>% .[N>=opts$min_gpc]

# Load RNA data
rna.sce <- load_SingleCellExperiment(io$rna.sce, normalise = TRUE, cells = opts$rna_cells)[args$rna_gene]

# Rename genes
# gene_metadata <- fread(io$gene.metadata) %>% .[ens_id%in%rownames(rna.sce) & symbol!=""]
# foo <- gene_metadata$symbol
# names(foo) <- gene_metadata$ens_id
# rna.sce <- rna.sce[rownames(rna.sce) %in% names(foo),]
# rownames(rna.sce) <- foo[rownames(rna.sce)]

# Extract data.table
rna.dt <- data.table(id_rna=colnames(rna.sce), value=logcounts(rna.sce)[1,], id=args$rna_gene)

# Merge data with sample metadata
acc.dt <- merge(acc.dt[,c("id_acc","id","value")], cell_metadata.dt[,c("id_acc","cell")], by="id_acc") %>% .[,id_acc:=NULL] %>% .[,assay:="Chromatin accessibility"]
met.dt <- merge(met.dt[,c("id_met","id","value")], cell_metadata.dt[,c("id_met","cell")], by="id_met") %>% .[,id_met:=NULL] %>% .[,assay:="DNA methylation"]
rna.dt <- merge(rna.dt[,c("id_rna","id","value")], cell_metadata.dt[,c("id_rna","cell")], by="id_rna") %>% .[,id_rna:=NULL] %>% .[,assay:="RNA expression"]

# bind in a single data table
to.plot <- do.call("rbind",list(rna.dt, met.dt, acc.dt)) %>% 
  merge(cell_metadata.dt[,c("cell","stage","celltype","class","celltype_class")], by="cell") %>% 
  .[,assay:=factor(assay,levels=c("RNA expression","DNA methylation","Chromatin accessibility"))]
  
##############
## Boxplots ##
##############

my_comparisons <- map(args$celltypes, function(i) sprintf('%s\n(%s)',i,c("WT","Tet-TKO")))

p <- ggplot(to.plot, aes(x=celltype_class, y=value)) +
  ggrastr::geom_jitter_rast(aes(fill=assay), shape=21, size=1.25, width=0.1, alpha=0.6, stroke=0.2) +
  geom_violin(aes(fill=assay), alpha=0.5, size=0.25) +
  geom_boxplot(aes(fill=assay), alpha=0.5, outlier.shape=NA, width=0.15, size=0.25) +
  stat_compare_means(comparisons = my_comparisons, method="t.test", size=2.5) +
  facet_wrap(~assay, ncol=1, scales="free_y") +
  # stat_summary(fun.data = give.n, geom = "text", size=2.5) +
  scale_fill_manual(values=opts$color) +
  scale_color_manual(values=opts$color) +
  labs(x="", y="", title="") +
  # guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    axis.title.y = element_text(colour="black", size=rel(1)),
    axis.text.x = element_text(size=rel(0.85), color="black"),
    axis.text.y = element_text(colour="black",size=rel(0.8)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position="none"
  )

pdf(sprintf("%s/plot_rna_%s_met_%s_acc_%s.pdf",args$outdir,args$rna_gene,sub(":","-",args$met_feature),sub(":","-",args$acc_feature)), width=6, height=5)
print(p)
dev.off()
