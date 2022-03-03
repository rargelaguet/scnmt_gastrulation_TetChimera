here::here("metacc/boxplots_feature_level/boxplots_feature_level.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

###################
## Load settings ##
###################

# I/O
io$metadata <- file.path(io$basedir,"results_new/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
io$outdir  <- file.path(io$basedir,"results_new/metacc/boxplots_feature_level/supp_fig"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min_observations <- 50
opts$min_cells <- 10
opts$min_marker_score <- 0.75
opts$annos <- c("prom_2000_2000","multiome_peaks")

# opts$celltypes <- c("Surface_ectoderm","Haematoendothelial_progenitors","Blood_progenitors","Endothelium","Pharyngeal_mesoderm","ExE_mesoderm")
opts$celltypes <- c("Surface_ectoderm","Pharyngeal_mesoderm","ExE_mesoderm")

opts$rename.celltypes <- c(
  "Erythroid1" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "early_Erythroid" = "Erythroid",
  "late_Erythroid" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak",
  "Mixed_mesoderm" = "Nascent_mesoderm",
  "Allantois" = "ExE_mesoderm"
)

# opts$markers.to.plot <- c("Surface_ectoderm","Haematoendothelial_progenitors","Blood_progenitors","Endothelium","Pharyngeal_mesoderm","ExE_mesoderm")
opts$markers.to.plot <- c("Surface_ectoderm","Pharyngeal_mesoderm","ExE_mesoderm")

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>% 
  .[,class:=ifelse(grepl("WT",class),"WT","TET-TKO")] %>% .[,class:=factor(class,levels=c("WT","TET-TKO"))] %>%
  .[celltype.mapped%in%opts$celltypes & (pass_metQC==TRUE | pass_accQC==TRUE)]

opts$met_cells <- sample_metadata[pass_metQC==TRUE,id_met]
opts$acc_cells <- sample_metadata[pass_accQC==TRUE,id_acc]


table(sample_metadata$celltype.mapped,sample_metadata$class)
# sample_metadata[celltype.mapped=="late_Erythroid",celltype.mapped:="Erythroid"]

tmp <- sample_metadata[,sum(pass_metQC,na.rm=T), by=c("celltype.mapped","class")]
#######################
## Load marker peaks ##
#######################

io$markers_file <- "/Users/argelagr/data/gastrulation_multiome_10x/results_new/atac/archR/differential/PeakMatrix/markers/marker_peaks_lenient.txt.gz"

marker_peaks.dt <- fread(io$markers_file) %>%
  .[celltype=="Erythroid2",celltype:="Erythroid"] %>%
  # .[,celltype:=stringr::str_replace_all(celltype,opts$rename.celltypes)] %>%
  .[celltype%in%opts$markers.to.plot] %>%
  .[score>=opts$min_marker_score] %>%
  .[,.(score=mean(score)), by=c("celltype","idx")]

table(marker_peaks.dt$celltype)

#######################
## Load met/acc data ##
#######################

metacc.dt <- opts$annos %>% map(function(i) {
  print(i)
  met.dt <- fread(file.path(io$basedir,sprintf("processed/met/feature_level/%s.tsv.gz",i))) %>% 
    setnames(c("id_met","id","anno","Nmet","Ntotal","rate")) %>% 
    .[id_met%in%opts$met_cells]
  
  acc.dt <- fread(file.path(io$basedir,sprintf("processed/acc/feature_level/%s.tsv.gz",i))) %>% 
    setnames(c("id_acc","id","anno","Nmet","Ntotal","rate")) %>% 
    .[id_acc%in%opts$acc_cells]
  
  if (i=="multiome_peaks") {
    met.dt <- met.dt[id%in%unique(marker_peaks.dt$idx)]
    acc.dt <- acc.dt[id%in%unique(marker_peaks.dt$idx)]
  }
  
  # Add common cell identifier
  met.dt <- merge(met.dt, sample_metadata[,c("cell","id_met")], by="id_met") %>% 
    .[,id_met:=NULL] %>% .[,context:="CG"]
  acc.dt <- merge(acc.dt, sample_metadata[,c("cell","id_acc")], by="id_acc") %>% 
    .[,id_acc:=NULL] %>% .[,context:="GC"]
  
  # Merge
  rbind(met.dt,acc.dt) %>% .[,anno:=factor(i)]
  
}) %>% rbindlist


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

# markers.to.plot <- c("Erythroid2","Surface_ectoderm")
# celltype.to.plot <- "late_Erythroid"
# anno.order <- c("prom_2000_2000","Erythroid2","Surface_ectoderm")

for (i in opts$celltypes) {
  
  to.plot <- rbind(
    metacc.dt[anno!="multiome_peaks"],
    metacc.dt[anno=="multiome_peaks"] %>%
      merge(marker_peaks.dt[,c("celltype","idx")] %>% setnames(c("celltype_marker","id")), by="id", allow.cartesian=TRUE) %>% .[,anno:=celltype_marker] %>% .[,celltype_marker:=NULL]
  ) %>% .[,.(Nmet=sum(Nmet), Ntotal=sum(Ntotal)), by=c("cell","anno","context")] %>%
    .[Ntotal>=opts$min_observations] %>%
    .[,rate:=100*(Nmet/Ntotal)] %>%
    merge(sample_metadata[celltype.mapped==i,c("cell","class","sample","celltype.mapped")], by="cell") %>%
    merge(global_rates.dt, by=c("cell","context")) %>%
    .[,rate_norm:=rate/global_rate] %>% .[,anno:=factor(anno,levels=c("prom_2000_2000",opts$markers.to.plot))]
  
  # for viz purposes
  to.plot[rate_norm>=3,rate_norm:=3]
  
  p.met <- ggplot(to.plot[context=="CG"], aes(x = class, y = rate_norm, fill=class)) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.35), size = 1.15, alpha = 0.50, shape=21) +
    facet_wrap(~anno, nrow=1, scales="fixed") +
    geom_hline(yintercept=1, linetype="dashed") +
    labs(x="", y="Methylation levels\n(relative to background)") +
    stat_compare_means(comparisons = list(c("WT", "TET-TKO")), aes(label = paste0("p = ", ..p.format..)), size=2, method="t.test") +
    # scale_fill_manual(values=opts$context.colors) +
    # coord_cartesian(ylim=c(0,100)) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size=rel(0.6)),
      legend.position = "none",
      legend.title = element_blank(),
      # axis.text = element_text(color="black"),
      # axis.text.x = element_blank(),
      axis.title.y = element_text(color="black", size=rel(0.80)),
      axis.text.x = element_text(color="black", size=rel(0.80)),
      axis.text.y = element_text(color="black", size=rel(0.80)),
      axis.ticks.x = element_blank()
    )
  
  p.acc <- ggplot(to.plot[context=="GC"], aes(x = class, y = rate_norm, fill=class)) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.35), size = 1.15, alpha = 0.50, shape=21) +
    facet_wrap(~anno, nrow=1, scales="fixed") +
    labs(x="", y="Chr. accessibility levels\n(relative to background)") +
    geom_hline(yintercept=1, linetype="dashed") +
    stat_compare_means(comparisons = list(c("WT", "TET-TKO")), aes(label = paste0("p = ", ..p.format..)), size=2, method="t.test") +
    # scale_fill_manual(values=opts$context.colors) +
    # coord_cartesian(ylim=c(8,60)) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = "none",
      legend.title = element_blank(),
      # axis.text = element_text(color="black"),
      # axis.text.x = element_blank(),
      axis.title.y = element_text(color="black", size=rel(0.80)),
      axis.text.x = element_text(color="black", size=rel(0.80)),
      axis.text.y = element_text(color="black", size=rel(0.80)),
      axis.ticks.x = element_blank()
    )
  
  pdf(file.path(io$outdir,sprintf("boxplots_metacc_fig_%s.pdf",i)), width=8, height=5)
  print(cowplot::plot_grid(plotlist=list(p.met,p.acc), nrow = 2))
  dev.off()
}


