here::i_am("rna/celltype_proportions/plot_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results_new/rna/celltype_proportions/fig")
dir.create(io$outdir, showWarnings = F)

# Options
opts$samples <- c(
  "E7.5_WT",
  "E7.5_TET_TKO",
  "E8.5_WT",
  "E8.5_WT_CD41+",
  "E8.5_WT_KDR+",
  "E8.5_WT_KDR+_CD41+",
  "E8.5_TET_TKO",
  "E8.5_TET_TKO_KDR+",
  "E8.5_TET_TKO_CD41+",
  "E8.5_TET_TKO_KDR+_CD41+"
)

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  # "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors",
  # "Blood_progenitors_2",
  # "Blood_progenitors_2",
  # "Erythroid1",
  # "Erythroid2",
  # "Erythroid3",
  "early_Erythroid",
  "late_Erythroid",
  "NMP",
  "Neurectoderm",
  # "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$rename <- c(
  # "Haematoendothelial_progenitors" = "Haemato._progenitors",
  "Caudal_epiblast" = "Epiblast"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & !is.na(celltype.mapped) & sample%in%opts$samples & celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename)] %>%
  setnames("celltype.mapped","celltype")

########################
## Barplot per sample ##
########################

# Calculate cell type proportions per sample
celltype_proportions.dt <- sample_metadata %>%
  .[,N:=.N,by="sample"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","celltype")] %>%
  setorder(sample)  %>% .[,sample:=factor(sample,levels=opts$samples)]

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(celltype_proportions.dt$celltype)]
celltype_proportions.dt[,celltype:=factor(celltype, levels=rev(names(opts$celltype.colors)))]

to.plot <- celltype_proportions.dt

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black") +
  scale_fill_manual(values=opts$celltype.colors) +
  facet_wrap(~sample, nrow=2, scales="free_x") +
  coord_flip() +
  labs(y="Number of cells") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    # strip.background =element_rect(fill=alpha("#DBDBDB", 0.50)),
    strip.text = element_text(color="black", size=rel(0.85)),
    axis.title.x = element_text(color="black", size=rel(0.9)),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=rel(1), color="black"),
    axis.text.x = element_text(size=rel(1), color="black")
  )

pdf(file.path(io$outdir,"celltype_proportions_fig.pdf"), width=12, height=7)
print(p)
dev.off()


#######################
## Barplot per class ##
#######################

to.plot <- sample_metadata %>%
  .[,class2:=ifelse(grepl("WT",class),"WT","Tet-TKO")] %>%
  .[,N:=.N,by="class2"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("class2","celltype")] %>%
  .[,class2:=factor(class2,levels=c("WT","Tet-TKO"))]

# Filter celltypes with small N
to.plot <- to.plot %>% .[,celltype:=factor(celltype,levels=rev(opts$celltypes[opts$celltypes%in%unique(celltype)]))]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black") +
  scale_fill_manual(values=opts$celltype.colors) +
  facet_wrap(~class2, nrow=1, scales="free_x") +
  coord_flip() +
  labs(y="Number of cells") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    # strip.background =element_rect(fill=alpha("#DBDBDB", 0.50)),
    strip.text = element_text(color="black", size=rel(0.85)),
    axis.title.x = element_text(color="black", size=rel(0.9)),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=rel(1), color="black"),
    axis.text.x = element_text(size=rel(1), color="black")
  )

pdf(file.path(io$outdir,"celltype_proportions_WT_vs_TKO.pdf"), width=5, height=6)
print(p)
dev.off()