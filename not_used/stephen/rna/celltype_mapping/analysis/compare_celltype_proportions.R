
#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

sample_metadata[, sample := plate]

# I/O
io$outdir <- paste0(io$basedir,"/results/celltype_proportions")
dir.create(paste0(io$outdir,"/boxplots"), showWarnings = F, recursive = TRUE)
dir.create(paste0(io$outdir,"/boxplots/per_class"), showWarnings = F)
dir.create(paste0(io$outdir,"/boxplots/per_sample"), showWarnings = F)

####################
## Define options ##
####################

opts$classes <- c(
 # "E7.5_chimera_KO",
#  "E7.5_crisprKO",
  "E8.5_host_WT",
  "E8.5_chimera_KO"
)

#opts$wt.classes <- opts$classes[opts$classes %like% "WT"]


  

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

opts$to.merge <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Intermediate_mesoderm" = "Mixed_mesoderm",
  "Paraxial_mesoderm" = "Mixed_mesoderm",
  "Nascent_mesoderm" = "Mixed_mesoderm",
  "Pharyngeal_mesoderm" = "Mixed_mesoderm"
  # "Visceral_endoderm" = "ExE_endoderm"
)


opts$remove.ExE.celltypes <- FALSE
opts$remove.blood <- FALSE
opts$remove.small.lineages <- FALSE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE ] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$to.merge)] %>% 
  .[, sample := plate]

names(opts$celltype.colors) <- stringr::str_replace_all(names(opts$celltype.colors), opts$to.merge)


# Filter cells
if (opts$remove.blood) {
  sample_metadata <- sample_metadata %>% .[!celltype.mapped=="Erythroid"]
}

if (opts$remove.ExE.celltypes) {
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

if (opts$remove.small.lineages) {
  opts$min.cells <- 25
  sample_metadata <- sample_metadata %>%
    .[,N:=.N,by=c("celltype.mapped")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
}

# Print statistics
table(sample_metadata$sample)
table(sample_metadata$celltype.mapped)


# just the raw numbers....
# plates <- sample_metadata[, .(plate, ko_type, stage, tdTOM, `KDR-Cy7`, `CD41-BV421`, class)] %>%
#   unique() %>% 
#   .[`KDR-Cy7` == TRUE, class := paste0(class, "_KDR+")] %>% 
#   .[`KDR-Cy7` == FALSE, class := paste0(class, "_KDR-")] %>%
#   .[`CD41-BV421` == TRUE, class := paste0(class, "_CD41+")] %>%
#   .[`CD41-BV421` == FALSE, class := paste0(class, "_CD41-")]

plates <- sample_metadata[plate %like% "E8.5", .(plate, ko_type, stage, tdTOM, `KDR-Cy7`, `CD41-BV421`, class)] %>%
  unique() %>% 
  .[, markers := paste0(`KDR-Cy7`, `CD41-BV421`)] %>% 
  .[markers == "FALSEFALSE", markers := "unsorted"] %>% 
  .[markers == "TRUETRUE", markers := "+/+"] %>% 
  .[markers == "FALSETRUE", markers := "CD41"] %>% 
  .[markers == "TRUEFALSE", markers := "_KDR"] %>%
  .[, ko := strsplit(class, "_") %>% map_chr(3)]

rawplot <- sample_metadata[, .N, .(plate, celltype.mapped)] %>%
  merge(plates[, .(plate, class, markers, ko)], by = "plate")

rawplot <- rawplot[class %like% "E8.5"]

sample_metadata[, table(class)]
plates[, table(class)]
rawplot[, table(class)]

sample_metadata[class %like% "E8.5", .N, .(class, `KDR-Cy7`, `CD41-BV421`, plate)][order(plate)]

p <- ggplot(rawplot, aes(x=celltype.mapped, y = N)) +
  geom_bar(aes(fill = celltype.mapped), stat="identity", alpha=0.9) +
  coord_flip()+#(ylim=c(-ylim,ylim)) +
  # geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  #scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  
  # labs(y="Difference in proportions (log2)", x="") +
  theme(
    legend.position = "none",
    # axis.title = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_fill_manual(
    breaks = names(opts$celltype.colors), 
    values = opts$celltype.colors
    ) +
  facet_wrap(ko~markers, nrow = 2)
p
####################################
## Calculate celltype proportions ##
####################################

# Calculate celltype proportions for WT
# wt_proportions.dt <- sample_metadata %>%
#   .[class%in%opts$wt.classes] %>%
#   .[,ncells:=.N, by="sample"] %>%
#   .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","class")]# %>%
# # .[,.(proportion=mean(proportion), N=mean(N)),by="celltype.mapped"]
wt_proportions.dt <- sample_metadata %>%
  .[class%in%opts$wt.classes] %>%
  .[,ncells:=.N] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped")]

# Calculate celltype proportions for KO samples
ko_proportions_per_sample.dt <- sample_metadata %>%
  .[!class%in%opts$wt.classes] %>%
  .[,ncells:=.N, by="sample"] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","sample","class")]

# Calculate celltype proportions for KO classes
ko_proportions_per_class.dt <- sample_metadata %>%
  .[!class%in%opts$wt.classes] %>%
  .[,ncells:=.N, by="class"] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","class")]

# Merge
proportions_per_sample.dt <- merge(
  ko_proportions_per_sample.dt, 
  wt_proportions.dt, 
  by = c("celltype.mapped"), allow.cartesian=T, suffixes = c(".ko",".wt")
) %>% .[,c("diff_proportion","diff_N"):=list(log2(proportion.ko/proportion.wt),N.ko-N.wt)] 

proportions_per_class.dt <- merge(
  ko_proportions_per_class.dt, 
  wt_proportions.dt, 
  by = c("celltype.mapped"), allow.cartesian=T, suffixes = c(".ko",".wt")
) %>% .[,c("diff_proportion","diff_N"):=list(log2(proportion.ko/proportion.wt),N.ko-N.wt)] 

#########################
## Barplots per sample ##
#########################

ylim <- max(abs(proportions_per_sample.dt$diff_proportion))

for (i in unique(proportions_per_sample.dt$sample)) {
  
  celltype.order <- proportions_per_sample.dt[sample==i,mean(diff_proportion),by="celltype.mapped"] %>% setorder(-V1) %>% .$celltype.mapped
  to.plot <- proportions_per_sample.dt[sample==i] %>% 
    .[,celltype.mapped:=factor(celltype.mapped,levels=celltype.order)]
  
  p <- ggplot(to.plot, aes(x=factor(celltype.mapped), y=diff_proportion)) +
    geom_point(aes(fill = celltype.mapped), shape=21, size=1) +
    geom_bar(aes(fill = celltype.mapped), stat="identity", alpha=0.5) +
    geom_text(y=-ylim, aes(label=N.wt), size=2.5) +
    geom_text(y=ylim, aes(label=N.ko), size=2.5) +
    coord_flip(ylim=c(-ylim,ylim)) +
    geom_hline(yintercept=0, linetype="dashed", size=0.5) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    theme_classic() +
    labs(y="Difference in proportions (log2)", x="") +
    theme(
      legend.position = "none",
      # axis.title = element_blank(),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black")
    )
  
  pdf(sprintf("%s/boxplots/per_sample/%s_barplots.pdf",io$outdir,i))
  print(p)
  dev.off()
}



# just the raw numbers....
plates <- sample_metadata[, .(plate, ko_type, stage, tdTOM, `KDR-Cy7`, `CD41-BV421`, class)] %>%
  unique() %>% 
  .[`KDR-Cy7` == TRUE, class := paste0(class, "_KDR+")] %>% 
  .[`KDR-Cy7` == FALSE, class := paste0(class, "_KDR-")] %>%
  .[`CD41-BV421` == TRUE, class := paste0(class, "_CD41+")] %>%
  .[`CD41-BV421` == FALSE, class := paste0(class, "_CD41-")]

rawplot <- sample_metadata[, .N, .(plate, celltype.mapped)] %>%
  merge(plates[, .(plate, class)], by = "plate")

p <- ggplot(rawplot, aes(x=celltype.mapped, y = N)) +
  geom_bar(aes(fill = celltype.mapped), stat="identity", alpha=0.5) +
  coord_flip()+#(ylim=c(-ylim,ylim)) +
  # geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  
  # labs(y="Difference in proportions (log2)", x="") +
  theme(
    legend.position = "none",
    # axis.title = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  facet_wrap(~class)


########################
## Boxplots per class ##  
#######################

ylim <- max(abs(proportions_per_sample.dt$diff_proportion))

for (i in unique(proportions_per_sample.dt$class)) {
  
  celltype.order <- proportions_per_sample.dt[class==i,mean(diff_proportion),by="celltype.mapped"] %>% setorder(-V1) %>% .$celltype.mapped
  to.plot <- proportions_per_sample.dt %>% 
    .[class==i] %>% 
    .[,celltype.mapped:=factor(celltype.mapped,levels=celltype.order)]
  
  tmp <- proportions_per_class.dt %>% 
    .[class==i] %>% 
    .[,celltype.mapped:=factor(celltype.mapped,levels=celltype.order)]
  
  p <- ggplot(to.plot, aes(x=factor(celltype.mapped), y=diff_proportion)) +
    geom_point(aes(fill = celltype.mapped), shape=21, size=1) +
    geom_boxplot(aes(fill = celltype.mapped), alpha=0.5) +
    geom_text(y=-ylim, aes(label=N.wt), size=2.5, data=tmp) +
    geom_text(y=ylim, aes(label=N.ko), size=2.5, data=tmp) +
    coord_flip(ylim=c(-ylim,ylim)) +
    geom_hline(yintercept=0, linetype="dashed", size=0.5) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    theme_classic() +
    labs(y="Difference in proportions (log2)", x="") +
    theme(
      legend.position = "none",
      # axis.title = element_blank(),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black")
    )
  
  pdf(sprintf("%s/boxplots/per_class/%s_boxplots.pdf",io$outdir,i))
  print(p)
  dev.off()
}

############################
## Polar plots per sample ##
############################

# %>% .[,diff_N_norm:=linMap(abs(diff_N), from=0.15, to=1.5)]

to.plot.wt_line <- data.table(
  celltype.mapped = unique(proportions_per_sample.dt$celltype.mapped),
  diff_proportion = log2(1),
  diff_N = 0 #+  0.01
)

p <- ggplot(proportions_per_sample.dt, aes(x=factor(celltype.mapped), y=diff_proportion, group=1)) +
  # geom_point(aes(color = celltype.mapped, size = diff_N_norm), stat = 'identity') + 
  geom_point(aes(color = celltype.mapped), size=2.5, stat = 'identity') + 
  geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line) +
  facet_wrap(~sample) +
  scale_color_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
  scale_size(guide = 'none') +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text =  element_text(size=rel(0.5)),
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank()
  )

pdf(sprintf("%s/polar_plots_by_sample.pdf",io$outdir), width=8, height=7)
print(p)
dev.off()



###########################
## Polar plots per class ##
###########################

# proportions_per_class.dt %>% .[,diff_N_norm:=linMap(abs(diff_N), from=0.15, to=1.5)]

to.plot.wt_line <- data.table(
  celltype.mapped = unique(proportions_per_class.dt$celltype.mapped),
  diff_proportion = log2(1),
  diff_N = 0 #+  0.01
)

p <- ggplot(proportions_per_class.dt, aes(x=factor(celltype.mapped), y=diff_proportion, group=1)) +
  # geom_point(aes(color = celltype.mapped, size = diff_N_norm), stat = 'identity') + 
  geom_point(aes(color = celltype.mapped), size=2.5, stat = 'identity') + 
  geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line) +
  facet_wrap(~class) +
  scale_color_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
  scale_size(guide = 'none') +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text =  element_text(size=rel(0.9)),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot.test$celltype.mapped)) * seq_along(to.plot.test$celltype.mapped))
  )

pdf(sprintf("%s/polar_plots_by_class.pdf",io$outdir), width=8, height=7)
print(p)
dev.off()

