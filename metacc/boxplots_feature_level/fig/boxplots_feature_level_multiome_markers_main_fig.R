here::here("metacc/boxplots_feature_level/boxplots_feature_level.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

###################
## Load settings ##
###################

# I/O
io$metadata <- file.path(io$basedir,"results_new/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
io$outdir  <- file.path(io$basedir,"results_new/metacc/boxplots_feature_level/fig")
dir.create(io$outdir, showWarnings = F)

# Options
opts$min_observations <- 50
opts$min_cells <- 10
opts$min_marker_score <- 0.75
opts$annos <- c("prom_2000_2000","LTR","CGI","multiome_peaks")

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>% 
  .[,class:=ifelse(grepl("WT",class),"WT","TET-TKO")] %>% .[,class:=factor(class,levels=c("WT","TET-TKO"))] %>%
  .[!is.na(celltype.mapped) & (pass_metQC==TRUE | pass_accQC==TRUE)]

opts$met_cells <- sample_metadata[pass_metQC==TRUE,id_met]
opts$acc_cells <- sample_metadata[pass_accQC==TRUE,id_acc]

###################
## Load TSS data ##
###################

metacc.dt <- opts$annos %>% map(function(i) {
  print(i)
  met.dt <- fread(file.path(io$basedir,sprintf("processed/met/feature_level/%s.tsv.gz",i))) %>% 
    setnames(c("id_met","id","anno","Nmet","Ntotal","rate")) %>% 
    .[id_met%in%opts$met_cells]
  
  acc.dt <- fread(file.path(io$basedir,sprintf("processed/acc/feature_level/%s.tsv.gz",i))) %>% 
    setnames(c("id_acc","id","anno","Nmet","Ntotal","rate")) %>% 
    .[id_acc%in%opts$acc_cells]
  
  # Add common cell identifier
  met.dt <- merge(met.dt, sample_metadata[,c("cell","id_met")], by="id_met") %>% 
    .[,id_met:=NULL] %>% .[,context:="CG"]
  acc.dt <- merge(acc.dt, sample_metadata[,c("cell","id_acc")], by="id_acc") %>% 
    .[,id_acc:=NULL] %>% .[,context:="GC"]
  
  # Merge
  metacc.dt <- rbind(met.dt,acc.dt) %>% .[,anno:=factor(i)]
  
}) %>% rbindlist


#######################
## Load marker peaks ##
#######################

io$markers_file <- "/Users/argelagr/data/gastrulation_multiome_10x/results_new/atac/archR/differential/PeakMatrix/markers/marker_peaks_lenient.txt.gz"
marker_peaks.dt <- fread(io$markers_file) 

scnmt.peaks <- unique(metacc.dt$id)
multiome.peaks <- unique(marker_peaks.dt$idx)

length(scnmt.peaks)
length(multiome.peaks)

# length(intersect(unique(marker_peaks.dt$idx),unique(metacc.dt$id)))

############################
## Calculate global rates ##
############################

global_rates.dt <- sample_metadata[,c("cell","met_rate","acc_rate")] %>%
  setnames(c("cell","CG","GC")) %>%
  melt(id.vars="cell", variable.name="context", value.name="global_rate")

# global_rates.dt <- metacc.dt %>% 
#   .[,.(Nmet=sum(Nmet), Ntotal=sum(Ntotal)), by=c("cell","anno","context")] %>% 
#   .[,global_rate:=100*(Nmet/Ntotal)] %>% .[,c("Nmet","Ntotal"):=NULL]

##########
## Plot ##
##########

markers.to.plot <- c("Erythroid2","Surface_ectoderm")
celltype.to.plot <- "late_Erythroid"

anno.order <- c("prom_2000_2000","CGI","LTR","Erythroid2","Surface_ectoderm")

to.plot <- rbind(
  metacc.dt[anno!="multiome_peaks"],
  metacc.dt[anno=="multiome_peaks"] %>%
    merge(marker_peaks.dt[celltype%in%markers.to.plot & score>=opts$min_marker_score,c("celltype","idx")] %>% setnames(c("celltype_marker","id")), by="id", allow.cartesian=TRUE) %>% .[,anno:=celltype_marker] %>% .[,celltype_marker:=NULL]
) %>% .[,.(Nmet=sum(Nmet), Ntotal=sum(Ntotal)), by=c("cell","anno","context")] %>%
  .[Ntotal>=opts$min_observations] %>%
  .[,rate:=100*(Nmet/Ntotal)] %>%
  merge(sample_metadata[celltype.mapped==celltype.to.plot,c("cell","class","sample","celltype.mapped")], by="cell") %>%
  merge(global_rates.dt, by=c("cell","context")) %>%
  .[,rate_norm:=rate/global_rate] %>%
  .[,anno:=factor(anno,levels=anno.order)]

p.met <- ggplot(to.plot[context=="CG"], aes(x = class, y = rate_norm, fill=class)) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35), size = 1.15, alpha = 0.50, shape=21) +
  facet_wrap(~anno, nrow=1, scales="fixed") +
  geom_hline(yintercept=1, linetype="dashed") +
  labs(x="", y="Methylation levels\n(relative to background)") +
  stat_compare_means(comparisons = list(c("WT", "TET-TKO")), aes(label = paste0("p = ", ..p.format..)), size=3, method="t.test") +
  scale_fill_manual(values=opts$class.colors) +
  # coord_cartesian(ylim=c(0,100)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    # axis.text = element_text(color="black"),
    # axis.text.x = element_blank(),
    axis.title.y = element_text(color="black", size=rel(0.80)),
    axis.text.x = element_text(color="black", size=rel(0.80)),
    axis.text.y = element_text(color="black", size=rel(0.80)),
    axis.ticks.x = element_blank()
  )

to.plot[context=="GC" & rate_norm>=2.5,rate_norm:=2.5] # for viz purposes
to.plot[context=="GC" & rate_norm<=0.70,rate_norm:=0.70] # for viz purposes

p.acc <- ggplot(to.plot[context=="GC"], aes(x = class, y = rate_norm, fill=class)) +
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35), size = 1.15, alpha = 0.50, shape=21) +
  facet_wrap(~anno, nrow=1, scales="fixed") +
  labs(x="", y="Chr. accessibility levels\n(relative to background)") +
  geom_hline(yintercept=1, linetype="dashed") +
  stat_compare_means(comparisons = list(c("WT", "TET-TKO")), aes(label = paste0("p = ", ..p.format..)), size=3, method="t.test") +
  scale_fill_manual(values=opts$class.colors) +
  # coord_cartesian(ylim=c(8,60)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    # axis.text = element_text(color="black"),
    # axis.text.x = element_blank(),
    axis.title.y = element_text(color="black", size=rel(0.80)),
    axis.text.x = element_text(color="black", size=rel(0.80)),
    axis.text.y = element_text(color="black", size=rel(0.80)),
    axis.ticks.x = element_blank()
  )

pdf(file.path(io$outdir,"boxplots_metacc_fig_legend.pdf"), width=6, height=5)
print(cowplot::plot_grid(plotlist=list(p.met,p.acc), nrow = 2))
dev.off()


