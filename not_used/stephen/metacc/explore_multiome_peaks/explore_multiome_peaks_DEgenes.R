## ----echo=TRUE, include=FALSE--------------------------------------------------------------------------------------------
source(here::here("settings.R"))


## ------------------------------------------------------------------------------------------------------------------------
io$outdir <- paste0(io$basedir,"/metacc/exploration/multiome_peaks")


## ------------------------------------------------------------------------------------------------------------------------


io$features.tsv <- io$features.tsv
io$DEgene_peaks <- gsub("\\..*", "_DEgenes.tsv.gz", io$features.tsv)

## ------------------------------------------------------------------------------------------------------------------------

# which genomic contexts?
opts$annos <- c(
  "multiome_peaks"
)

# which DE gene comparisons?
opts$DEcomparisons <- c(
  "Erythroid1_E8.5_Host_vs_E8.5_TET_TKO",
  "Erythroid2_E8.5_Host_vs_E8.5_TET_TKO",
  "Erythroid3_E8.5_Host_vs_E8.5_TET_TKO"
)
## ------------------------------------------------------------------------------------------------------------------------

celltypes <- sample_metadata[, unique(celltype.mapped)]


## ------------------------------------------------------------------------------------------------------------------------
acc.stats <- fread(io$acc.stats) %>% .[,c("id_acc","mean")] %>% .[,context:="GC"]
met.stats <- fread(io$met.stats) %>% .[,c("id_met","mean")] %>% .[,context:="CG"]


## ------------------------------------------------------------------------------------------------------------------------

marker_peaks.dt <- fread(io$DEgene_peaks) %>% 
  .[celltype %in% opts$DEcomparisons] %>% 
  .[, c("idx", "id") := .(id, sprintf("%s:%s-%s",chr,start,end))]



head(marker_peaks.dt)



## ------------------------------------------------------------------------------------------------------------------------

feature_metadata <- fread(io$features.tsv) %>% 
  setnames(c("chr", "start", "end", "strand", "id", "anno"))

# correct chromosome naming
if (!feature_metadata[1, chr %like% "chr"]) feature_metadata[, chr := paste0("chr", chr)]

feature_metadata[,idx:=sprintf("%s:%s-%s",chr,start,end)]

head(feature_metadata)


## ----load_metdata, echo=FALSE, include=FALSE-----------------------------------------------------------------------------

met_dt <- paste0(io$met_data_parsed, "/", opts$annos, ".tsv.gz") %>% 
  map(fread) %>% 
  rbindlist()

head(met_dt)


met_metadata <- sample_metadata[pass_metQC == TRUE] %>% 
  .[, cell := gsub(".tsv.gz", "", basename(cg_files))]

to.plot_met_all <- merge(
  met_dt,
  met_metadata,
  by = "cell"
) %>% 
  .[, .(rate = mean(rate)), .(cell, stage_lineage, anno, class, markers )]


## ----load_accdata, echo=FALSE, include=FALSE-----------------------------------------------------------------------------

acc_dt <- paste0(io$acc_data_parsed, "/", opts$annos, ".tsv.gz") %>% 
  map(fread) %>% 
  rbindlist()

head(acc_dt)





## ------------------------------------------------------------------------------------------------------------------------
met_dt.filt <- merge(
  met_dt,
  marker_peaks.dt[,c("id","celltype","updown")], 
  by="id"
)

head(met_dt.filt)

acc_dt.filt <- merge(
  acc_dt, 
  marker_peaks.dt[,c("id","celltype","updown")], 
  by="id"
)

head(met_dt.filt)

## ------------------------------------------------------------------------------------------------------------------------
length(unique(met_dt.filt$id))
length(unique(acc_dt.filt$id))


## ------------------------------------------------------------------------------------------------------------------------

met_metadata <- sample_metadata[pass_metQC == TRUE] %>% 
  .[, cell := gsub(".tsv.gz", "", basename(cg_files))]

to.plot_met <- merge(
  met_dt.filt,
  met_metadata,
  by = "cell"
) %>% 
  .[, .(rate = mean(rate)), .(cell, stage_lineage, anno, celltype, class, markers, updown)]




to.plot_met <- merge(
  to.plot_met, 
  met.stats[, .(cell = id_met, mean)],
  by = "cell"
) %>% 
  .[, scaled := rate / (100*mean)]

##########################################################

# repeat for accessibility

acc_metadata <- sample_metadata[pass_accQC == TRUE] %>% 
  .[, cell := gsub(".tsv.gz", "", basename(gc_files))]

to.plot_acc <- merge(
  acc_dt.filt,
  acc_metadata,
  by = "cell"
) %>% 
  .[, .(rate = mean(rate)), .(cell, stage_lineage, anno, celltype, class, markers)]




to.plot_acc <- merge(
  to.plot_acc, 
  acc.stats[, .(cell = id_acc, mean)],
  by = "cell"
) %>% 
  .[, scaled := rate / (100*mean)]







## ------------------------------------------------------------------------------------------------------------------------


# features_sub <- c(
#   "Cardiomyocytes",
#   "Erythroid",
#   "Forebrain_Midbrain_Hindbrain" ,
#   "Gut",
#   "ExE_endoderm"
# )

lineages_sub <- to.plot_met[, length(unique(cell)), stage_lineage] %>% 
  .[V1 >= 20, stage_lineage]

to.plot_met_all[, celltype := "All Peaks"]

p <- ggplot(to.plot_met_all[stage_lineage %in% lineages_sub], aes(x = stage_lineage, y = rate, fill = class)) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.1, alpha = 0.2) +
  labs(x="", y="Methylation") +
  facet_wrap(~celltype, scales="fixed") +
  #scale_fill_manual(values=opts$colors) +
  # coord_cartesian(ylim=c(5,60)) +
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    #legend.position = "none",
    axis.text = element_text(color="black")
  )
p


p <- ggplot(to.plot_met[stage_lineage %in% lineages_sub], aes(x = stage_lineage, y = rate, fill = class)) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.1, alpha = 0.2) +
  labs(x="", y="Methylation") +
  facet_wrap(updown~celltype, scales="fixed") +
  #scale_fill_manual(values=opts$colors) +
  # coord_cartesian(ylim=c(5,60)) +
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    #legend.position = "none",
    axis.text = element_text(color="black")
  )
p

dir.create(io$outdir, recursive = TRUE)
cowplot::save_plot(paste0(io$outdir, "/met_boxplots_celltypes.pdf"), p, base_height = 20, base_width = 20)


p <- ggplot(to.plot_met[celltype %in% features_sub & stage_lineage %in% lineages_sub], aes(x = markers, y = rate, fill = class)) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.1, alpha = 0.2) +
  labs(x="", y="Methylation") +
  facet_wrap(~celltype, scales="fixed") +
  #scale_fill_manual(values=opts$colors) +
  # coord_cartesian(ylim=c(5,60)) +
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    #legend.position = "none",
    axis.text = element_text(color="black")
  )
p

dir.create(io$outdir, recursive = TRUE)
cowplot::save_plot(paste0(io$outdir, "/met_boxplots_markers.pdf"), p, base_height = 20, base_width = 20)




## ------------------------------------------------------------------------------------------------------------------------


p <- ggplot(to.plot_acc[celltype %in% features_sub & stage_lineage %in% lineages_sub], aes(x = stage_lineage, y = rate, fill = class)) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.1, alpha = 0.2) +
  labs(x="", y="Accessibility") +
  facet_wrap(~celltype, scales="fixed") +
  #scale_fill_manual(values=opts$colors) +
  # coord_cartesian(ylim=c(5,60)) +
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    #legend.position = "none",
    axis.text = element_text(color="black")
  )
p

dir.create(io$outdir, recursive = TRUE)
cowplot::save_plot(paste0(io$outdir, "/acc_boxplots_celltypes.pdf"), p, base_height = 20, base_width = 20)


p <- ggplot(to.plot_acc[celltype %in% features_sub & stage_lineage %in% lineages_sub], aes(x = markers, y = rate, fill = class)) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.1, alpha = 0.2) +
  labs(x="", y="Accessibility") +
  facet_wrap(~celltype, scales="fixed") +
  #scale_fill_manual(values=opts$colors) +
  # coord_cartesian(ylim=c(5,60)) +
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    #legend.position = "none",
    axis.text = element_text(color="black")
  )
p

dir.create(io$outdir, recursive = TRUE)
cowplot::save_plot(paste0(io$outdir, "/acc_boxplots_markers.pdf"), p, base_height = 20, base_width = 20)




## ------------------------------------------------------------------------------------------------------------------------
# to.plot2 <- to.plot %>% merge(sample_metadata[,c("id_rna","lineage10x")])
# p <- ggplot(to.plot2[stage_lineage=="E7.5_Endoderm"], aes(x=lineage10x, y=rate, fill=lineage10x)) +
#   geom_boxplot(outlier.shape=NA, coef=1) +
#   labs(x="", y="Methylation (%)") +
#   facet_wrap(~sign, nrow=1, scales="fixed") +
#   # coord_cartesian(ylim=c(5,60)) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   theme_bw() +
#   guides(x = guide_axis(angle = 90)) +
#   theme(
#     legend.position = "none",
#     axis.text = element_text(color="black")
#   )
# print(p)


## ------------------------------------------------------------------------------------------------------------------------
acc.stats <- acc_dt[,.(mean=mean(rate)),by="id_acc"] %>% .[,context:="GC"]
met.stats <- met_dt[,.(mean=mean(rate)),by="id_met"] %>% .[,context:="CG"]


## ------------------------------------------------------------------------------------------------------------------------
theme_pub <- function() {
    theme(
      axis.text.x = element_text(size=rel(0.50), color="black"),
      axis.text.y = element_text(size=rel(0.75), color="black"),
      axis.title = element_text(size=rel(0.75), color="black"),
      legend.position = "none"
    )
}


## ------------------------------------------------------------------------------------------------------------------------
for (i in unique(marker_peaks.dt$celltype)) {
  
  ## Prepare data
  
  to.plot.met <- met_dt %>% 
    .[id%in%marker_peaks.dt[celltype==i,id]] %>%
    .[,.(rate=100*(sum(Nmet)/sum(Ntotal))),by=c("id_met")] %>%
    merge(sample_metadata[,c("id_met","id_rna","lineage10x","stage_lineage")], by="id_met") %>%
    merge(met.stats,by="id_met") %>% .[,norm_rate:=rate/mean]  %>%
    # .[,lineage10x:=factor(lineage10x, levels=opts$)]
    .[,stage_lineage:=factor(stage_lineage, levels=opts$stage_lineage)]
  
  to.plot.acc <- acc_dt %>% 
    .[id%in%marker_peaks.dt[celltype==i,id]] %>%
    .[,.(rate=100*(sum(Nmet)/sum(Ntotal))),by=c("id_acc")] %>%
    merge(sample_metadata[,c("id_acc","id_rna","lineage10x","stage_lineage")], by="id_acc") %>%
    merge(acc.stats,by="id_acc") %>% .[,norm_rate:=rate/mean]  %>%
    # .[,lineage10x:=factor(lineage10x, levels=opts$)]
    .[,stage_lineage:=factor(stage_lineage, levels=opts$stage_lineage)]
  
  ## Plot
  
  p.met.unscaled <- ggboxplot(to.plot.met, x = "stage_lineage", y = "rate", fill="stage_lineage", coef=0.5, outlier.shape=NA) +
    scale_fill_manual(values=opts$stagelineage.colors) +
    labs(x="", y="Methylation (%)") +
    guides(x = guide_axis(angle = 90)) +
    coord_cartesian(ylim=c(0,100)) +
    theme_pub()
  p.met.scaled <- ggboxplot(to.plot.met, x = "stage_lineage", y = "norm_rate", fill="stage_lineage", coef=0.5, outlier.shape=NA) +
    geom_hline(yintercept=1, linetype="dashed", alpha=0.75, size=0.75) +
    scale_fill_manual(values=opts$stagelineage.colors) +
    labs(x="", y="Methylation (% relative to bgd)") +
    guides(x = guide_axis(angle = 90)) +
    coord_cartesian(ylim=c(0.1,1.75)) +
    theme_pub()
  
  p.acc.unscaled <- ggboxplot(to.plot.acc, x = "stage_lineage", y = "rate", fill="stage_lineage", coef=0.5, outlier.shape=NA) +
    scale_fill_manual(values=opts$stagelineage.colors) +
    labs(x="", y="Accessibility (%)") +
    coord_cartesian(ylim=c(20,70)) +
    guides(x = guide_axis(angle = 90)) +
    theme_pub()
  
  p.acc.scaled <- ggboxplot(to.plot.acc, x = "stage_lineage", y = "norm_rate", fill="stage_lineage", coef=0.5, outlier.shape=NA) +
    geom_hline(yintercept=1, linetype="dashed", alpha=0.75, size=0.75) +
    scale_fill_manual(values=opts$stagelineage.colors) +
    labs(x="", y="Accessibility (% relative to bgd)") +
    coord_cartesian(ylim=c(0.6,1.8)) +
    guides(x = guide_axis(angle = 90)) +
    theme_pub()
  
  p <- cowplot::plot_grid(plotlist=list(p.met.unscaled, p.met.scaled, p.acc.unscaled, p.acc.scaled), nrow = 2, ncol=2)
  pdf(sprintf("%s/%s_markerPeaks_scnmt_metacc.pdf",io$outdir,i), width=11, height=8)
  print(p)
  dev.off()
}

