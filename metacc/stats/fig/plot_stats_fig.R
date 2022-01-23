here::i_am("metacc/stats/plot_stats.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Sanity checks
stopifnot(opts$context %in% c("CG","GC"))

# I/O
io$outdir <- file.path(io$basedir,"results_new/metacc/stats/pdf")
dir.create(io$outdir, showWarnings=F)

# Options
opts$context <- "CG"

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

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[celltype.mapped%in%opts$celltypes] %>% .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)]

# Merge with cell metadata
if (opts$context=="CG") {
  to.plot <- sample_metadata %>%
    .[!is.na(id_met)] %>%
    .[,c("cell","plate","sample","celltype.mapped","nCG","met_rate")] %>% 
    setnames(c("nCG","met_rate"),c("N","rate"))
} else if (opts$context=="GC") {
  to.plot <- sample_metadata %>%
    .[!is.na(id_acc)] %>%
    .[,c("cell","plate","sample","celltype.mapped","nGC","acc_rate")] %>% 
    setnames(c("nGC","acc_rate"),c("N","rate"))
}

######################################
## Boxplots with rate per cell type ##
######################################

to.plot2 <- to.plot %>% 
  .[,ko:=ifelse(grepl("KO",sample),"TET TKO","WT")] %>%
  .[!is.na(celltype.mapped)] %>% .[,N:=.N,by="celltype.mapped"] %>% .[N>=5] %>% 
  .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes[opts$celltypes%in%unique(celltype.mapped)])]

p <- ggboxplot(to.plot2, x = "celltype.mapped", y = "rate", outlier.shape=NA, fill="celltype.mapped", alpha=0.75) +
  geom_jitter(alpha=0.50, size=0.75, shape=21, width=0.1) +
  facet_wrap(~ko, nrow=2, scales="fixed") + 
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Global methylation levels (%)") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size=rel(0.70)),
    axis.title.y = element_text(size=rel(0.90)),
    axis.text.x = element_text(size=rel(0.65))
  )

pdf(file.path(io$outdir,sprintf("%s_global_levels.pdf",opts$context)), width=6, height=5)
print(p)
dev.off()

##################################
## Boxplots with rate per class ##
##################################

opts$celltypes = c(
  "Spinal_cord",
  # "Surface_ectoderm",
  "Gut",
  "Pharyngeal_mesoderm",
  "Endothelium",
  "Haematoendothelial_progenitors",
  "Blood_progenitors",
  "early_Erythroid",
  "late_Erythroid"
)

to.plot3 <- to.plot2[celltype.mapped%in%opts$celltypes]

p <- ggboxplot(to.plot3, x = "ko", y = "rate", outlier.shape=NA, fill="ko", alpha=0.75) +
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)), label.y = 89, method="t.test", size=3) +
  geom_jitter(alpha=0.50, size=0.75, shape=21, width=0.1) +
  facet_wrap(~celltype.mapped, nrow=2, scales="fixed") + 
  coord_cartesian(ylim=c(50,91)) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Global methylation levels (%)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size=rel(0.80), color="black"),
    axis.title.y = element_text(size=rel(0.90), color="black"),
    axis.text.x = element_text(size=rel(0.80), color="black")
  )

pdf(file.path(io$outdir,sprintf("%s_global_levels_wt_vs_ko.pdf",opts$context)), width=7, height=4.5)
print(p)
dev.off()


###################################################
## Number of cells with multi-omics measurements ##
###################################################

to.plot <- sample_metadata %>%
  .[,ko:=ifelse(grepl("KO",class),"Tet-TKO","WT")] %>%
  .[,.(met=sum(!is.na(id_met)), acc=sum(!is.na(id_acc)), rna=sum(!is.na(id_rna))),by="ko"] %>%
  melt(id.vars="ko")

p <- ggplot(to.plot, aes_string(x="ko", y="value", fill="variable")) +
  geom_bar(stat="identity", position="dodge", color="black") +
  labs(x="", y="Number of cells") +
  scale_fill_manual(values=c(acc="#00BFC4", rna="#00CD00", met="#F8766D")) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.text.y = element_text(size=rel(1.0), color="black"),
    axis.text.x = element_text(size=rel(1.2), color="black")
  )

pdf(file.path(io$outdir,"number_cells_per_omic.pdf"), width=7, height=5)
print(p)
dev.off()
