## ----echo=TRUE, include=FALSE--------------------------------------------------------------------------------------------
source(here::here("settings.R"))


## ------------------------------------------------------------------------------------------------------------------------
io$outdir <- paste0(io$basedir,"/metacc/exploration/global_boxplots")


## ------------------------------------------------------------------------------------------------------------------------

dir(io$features.dir)
# Define genomic contexts
opts$annos <- c(
  "100000bp_bins_50000_step"
  #"10000bp_bins_5000_step"
)

## ------------------------------------------------------------------------------------------------------------------------




celltypes <- sample_metadata[, unique(celltype.mapped)]


## ------------------------------------------------------------------------------------------------------------------------
acc.stats <- fread(io$acc.stats) %>% .[,c("id_acc","mean")] %>% .[,context:="GC"]
met.stats <- fread(io$met.stats) %>% .[,c("id_met","mean")] %>% .[,context:="CG"]



## ------------------------------------------------------------------------------------------------------------------------

feature_metadata <- paste0(io$features.dir, "/", opts$annos, ".bed.gz") %>% 
  map(fread) %>% 
  rbindlist() %>% 
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


## ----load_accdata, echo=FALSE, include=FALSE-----------------------------------------------------------------------------

acc_dt <- paste0(io$acc_data_parsed, "/", opts$annos, ".tsv.gz") %>% 
  map(fread) %>% 
  rbindlist()

head(acc_dt)







## ------------------------------------------------------------------------------------------------------------------------

met_metadata <- sample_metadata[pass_metQC == TRUE] %>% 
  .[, cell := gsub(".tsv.gz", "", basename(cg_files))]

to.plot_met <- merge(
  met_dt,
  met_metadata,
  by = "cell"
) %>% 
  .[, .(rate = mean(rate)), .(cell, stage_lineage, anno, class, markers)]

head(to.plot_met)



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
  acc_dt,
  acc_metadata,
  by = "cell"
) %>% 
  .[, .(rate = mean(rate)), .(cell, stage_lineage, anno,  class, markers)]




to.plot_acc <- merge(
  to.plot_acc, 
  acc.stats[, .(cell = id_acc, mean)],
  by = "cell"
) %>% 
  .[, scaled := rate / (100*mean)]







## ------------------------------------------------------------------------------------------------------------------------



lineages_sub <- to.plot_met[, length(unique(cell)), stage_lineage] %>% 
  .[V1 >= 10 & !stage_lineage %like% "NA", stage_lineage]


p <- ggplot(to.plot_met[stage_lineage %in% lineages_sub], aes(x = stage_lineage, y = rate, fill = class)) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.1, alpha = 0.2) +
  labs(x="", y="Methylation") +
  facet_wrap(~anno, scales="fixed") +
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


p <- ggplot(to.plot_met[ stage_lineage %in% lineages_sub], aes(x = markers, y = rate, fill = class)) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.1, alpha = 0.2) +
  labs(x="", y="Methylation") +
  facet_wrap(~anno, scales="fixed") +
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


p <- ggplot(to.plot_acc[ stage_lineage %in% lineages_sub], aes(x = stage_lineage, y = rate, fill = class)) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.1, alpha = 0.2) +
  labs(x="", y="Accessibility") +
  facet_wrap(~anno, scales="fixed") +
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


p <- ggplot(to.plot_acc[stage_lineage %in% lineages_sub], aes(x = markers, y = rate, fill = class)) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.1, alpha = 0.2) +
  labs(x="", y="Accessibility") +
  facet_wrap(~anno, scales="fixed") +
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


