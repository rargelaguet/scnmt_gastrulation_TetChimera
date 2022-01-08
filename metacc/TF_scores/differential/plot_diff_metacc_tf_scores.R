# here::i_am("")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results_new/metacc/TF_scores/differentialpdf"); dir.create(args$outdir, showWarnings=F)

# Options

# Minimum differential (%) for statistical significance
opts$min.diff <- 25

# Multiple testing correction
opts$threshold_fdr <- 0.10

###############
## Load data ##
###############

diff.dt <- fread("/Users/argelagr/data/tet_chimera_nmtseq/results_new/metacc/TF_scores/differential/KO_vs_WT.txt.gz")

#############################
## Plot diffmet vs diffacc ##
#############################

to.plot <- diff.dt %>% dcast(tf+celltype~context, value.var=c("diff","padj_fdr")) %>%
  .[complete.cases(.)]

to.plot[,dot_size:=minmax.normalisation(abs(diff_CG))]
to.plot.text <- to.plot %>% split(.$celltype) %>% map(~ sort.abs(.,"diff_CG") %>% head(n=10)) %>% rbindlist

ggscatter(to.plot, x="diff_CG", y="diff_GC", shape=21, size="dot_size", fill="dot_size") +
  facet_wrap(~celltype, scales="fixed") +
  labs(x="Methylation difference (%)", y="Accessibility difference (%)") +
  scale_size_continuous(range = c(0.25,2)) +
  ggrepel::geom_text_repel(data=to.plot.text, aes(x=diff_CG, y=diff_GC, label=tf), size=3,  max.overlaps=100, segment.color = NA) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_fill_gradient2(low = "gray50", mid="gray90", high = "red") +
  # coord_cartesian(xlim=c(0,0.9)) +
  theme_classic() +
  theme(
    axis.text = element_text(color="black"),
    # legend.position = "right",
    legend.position = "none",
    legend.title = element_blank()
  )

# pdf(file.path(io$outdir,"lineplot_number_binding_sites_per_tf_positive_vs_negative.pdf"), width=6, height=5)
# print(p)
# dev.off()


#############################
## Plot diffmet vs diffRNA ##
#############################

####################################
## Load pseudobulk RNA expression ##
####################################

sce <- readRDS("/Users/argelagr/data/tet_chimera_nmtseq/results_new/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype_ko.rds")

sce <- sce[,grep("WT",colnames(sce))]
colnames(sce) <- gsub("_WT","",colnames(sce))

TFs <- intersect(toupper(rownames(sce)),unique(diff.dt$tf))
tf.sce <- sce[str_to_title(TFs)]; rownames(tf.sce) <- toupper(rownames(tf.sce))

opts$rename.celltypes <- c(
  "Erythroid1" = "early_Erythroid",
  "Erythroid2" = "early_Erythroid",
  "Erythroid3" = "late_Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak",
  "Mixed_mesoderm" = "Nascent_mesoderm",
  "Allantois" = "ExE_mesoderm"
)

rna.dt <- as.matrix(logcounts(tf.sce)) %>% t %>%
  as.data.table(keep.rownames="celltype") %>% 
  melt(id.vars="celltype", value.name="expr", variable.name="tf") %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename.celltypes)] %>%
  .[,.(expr=mean(expr)),by=c("tf","celltype")]


####################################
## Plot diffmet vs RNA expression ##
####################################

to.plot <- diff.dt %>% 
  dcast(tf+celltype~context, value.var=c("diff","padj_fdr")) %>%
  .[complete.cases(.)] %>%
  merge(rna.dt,by=c("celltype","tf"))

to.plot[,dot_size:=minmax.normalisation(abs(diff_CG)*expr)]
to.plot.text <- to.plot %>% split(.$celltype) %>% map(~ sort.abs(.,"dot_size") %>% head(n=15)) %>% rbindlist

ggscatter(to.plot, x="diff_CG", y="expr", shape=21, size="dot_size", fill="dot_size") +
  facet_wrap(~celltype, scales="fixed") +
  labs(x="Methylation difference (%)", y="RNA expression") +
  scale_size_continuous(range = c(0.25,2.5)) +
  ggrepel::geom_text_repel(data=to.plot.text, aes(x=diff_CG, y=expr, label=tf), size=3,  max.overlaps=100, segment.color = NA) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_fill_gradient2(low = "gray50", mid="gray90", high = "red") +
  # coord_cartesian(xlim=c(0,0.9)) +
  theme_classic() +
  theme(
    axis.text = element_text(color="black"),
    # legend.position = "right",
    legend.position = "none",
    legend.title = element_blank()
  )

# pdf(file.path(io$outdir,"lineplot_number_binding_sites_per_tf_positive_vs_negative.pdf"), width=6, height=5)
# print(p)
# dev.off()
