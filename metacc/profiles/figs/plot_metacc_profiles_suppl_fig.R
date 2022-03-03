here::here("metacc/profiles/calculate_metacc_profiles.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",  help='Cell metadata')
# p$add_argument('--anno',    type="character",    help='Genomic annotation')
p$add_argument('--file',    type="character",    help='Precomputed file')
p$add_argument('--outdir',  type="character",    help='Output directory')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

# I/O
io$metadata <- file.path(io$basedir,"results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
io$outdir <- file.path(io$basedir,"results/metacc/profiles/supp_fig"); dir.create(io$outdir, showWarnings = F)

# Options
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
  .[,celltype:=stringr::str_replace_all(celltype.mapped,opts$rename.celltypes)] %>%
  .[,class:=ifelse(grepl("WT",class),"WT","TET-TKO")] %>% .[,class:=factor(class,levels=c("WT","TET-TKO"))] %>%
  .[(pass_metQC==TRUE | pass_accQC==TRUE) & celltype%in%opts$celltypes]

table(sample_metadata$celltype,sample_metadata$class)

# tmp <- fread(io$metadata) %>% .[,celltype:=stringr::str_replace_all(celltype.mapped,opts$rename.celltypes)] %>% .[,class:=ifelse(grepl("WT",class),"WT","TET-TKO")] %>% .[,class:=factor(class,levels=c("WT","TET-TKO"))] %>%.[(pass_metQC==TRUE | pass_accQC==TRUE)]
# table(tmp$celltype,tmp$class)

#######################
## Load marker peaks ##
#######################

opts$min_marker_score <- 0.75

marker_peaks.dt <- fread(io$multiome.marker_peaks) %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename.celltypes)] %>%
  .[celltype%in%opts$markers.to.plot] %>%
  # .[celltype%in%opts$celltypes] %>%
  .[score>=opts$min_marker_score] %>%
  .[,.(score=mean(score)), by=c("celltype","idx")] %>%
  .[,idx:=gsub("chr","",idx)]

table(marker_peaks.dt$celltype)

#######################
## Load TSS profiles ##
#######################

io$precomputed_metacc_profiles  <- file.path(io$basedir,"results_new/metacc/profiles/tss/precomputed_metacc_prom_200_200.txt.gz")
metacc_promoters.dt <- fread(io$precomputed_metacc_profiles) %>% .[cell%in%sample_metadata$cell]

stopifnot(sample_metadata$cell%in%unique(metacc_promoters.dt$cell))
sample_metadata[!cell%in%unique(metacc_promoters.dt$cell),c("cell","id_met","id_acc")]

#################################
## Load Multiome peak profiles ##
#################################

io$precomputed_metacc_profiles  <- file.path(io$basedir,"results_new/metacc/profiles/multiome_peaks/precomputed_metacc_multiome_peaks_filt_v2.txt.gz")
# io$precomputed_metacc_profiles  <- file.path(io$basedir,"results_new/metacc/profiles/multiome_peaks/first_trial/precomputed_metacc_multiome_peaks_filt.txt.gz")
metacc_multiome.dt <- fread(io$precomputed_metacc_profiles) %>%
  .[cell%in%sample_metadata$cell & id%in%unique(marker_peaks.dt$idx)]

# Sanity checks
# metacc_multiome.dt[,.N,by=c("cell","context")] %>% View
# stopifnot(sample_metadata$cell%in%unique(metacc_multiome.dt$cell))
# sample_metadata[!cell%in%unique(metacc_multiome.dt$cell),c("cell","id_met","id_acc")]

metacc_multiome.dt <- metacc_multiome.dt %>%
  merge(marker_peaks.dt[,c("celltype","idx")] %>% setnames(c("celltype_marker","id")), by="id", allow.cartesian=T) %>%
  .[,anno:=NULL] %>% setnames("celltype_marker","anno")

#############
## Combine ##
#############

metacc.dt <- rbind(metacc_promoters.dt,metacc_multiome.dt)

# metacc.dt[,c("id","anno")] %>% unique %>% .[,.N,by="anno"]
# rm(metacc_promoters.dt,metacc_multiome.dt)

###########################################
## Plot TSS profiles one class at a time ##
###########################################

# i <- "ExE_mesoderm"

for (i in opts$celltypes) {
  
  # Methylation
  to.plot.met <- metacc.dt %>%
    .[cell%in%sample_metadata[pass_metQC==TRUE & celltype==i,cell]] %>%
    .[,.(rate=mean(rate), N=sum(N)),by=c("dist","context","cell","anno")]
  
  to.plot.acc <- metacc.dt %>%
    .[cell%in%sample_metadata[pass_accQC==TRUE & celltype==i,cell]] %>%
    .[,.(rate=mean(rate), N=sum(N)),by=c("dist","context","cell","anno")]
  
  to.plot <- rbind(to.plot.met,to.plot.acc) %>% 
    merge(sample_metadata[,c("cell","class")]) %>%
    .[,anno:=factor(anno,levels=c("prom_200_200",opts$markers.to.plot))]
    
  # fILTER OUT CELLS WITH SMALL N
  to.plot <- to.plot[N>=15]
  
  # anno.order <- c("prom_200_200","Erythroid","Surface_ectoderm")
  
  to.plot.lines <- sample_metadata[celltype==i,.(CG=mean(met_rate,na.rm=T), GC=mean(acc_rate,na.rm=T)), by="class"] %>%
    melt(id.vars=c("class"), variable.name="context", value.name="rate")
  
  p <- ggplot(to.plot, aes(x=dist, y=rate, group=context, fill=context, color=context)) +
    stat_summary(geom="ribbon", fun.data="mean_sd", alpha=1, color="black") +
    # stat_summary(geom="line", fun.data="mean_se", size=0.75) +
    # geom_line(size=2) +
    facet_wrap(~class~anno, scales="fixed", nrow=2) +
    geom_hline(aes(yintercept=rate, color=context), linetype="dashed", alpha=0.75, size=0.75, data=to.plot.lines) +
    labs(x="Distance from center (bp)", y="Met/Acc levels (%)") +
    coord_cartesian(ylim=c(12,90)) +
    # scale_x_continuous(breaks=c(-1,0,1)) +
    # xlim(-opts$window_size, opts$window_size) +
    guides(fill="none", color="none", linetype="none") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size=rel(1.0), colour="black"),
      axis.text.y = element_text(size=rel(1.1), colour="black")
    )
  
  pdf(file.path(io$outdir,sprintf("metacc_profiles_%s_fig.pdf",i)), width=8, height=5)
  print(p)
  dev.off()
}

