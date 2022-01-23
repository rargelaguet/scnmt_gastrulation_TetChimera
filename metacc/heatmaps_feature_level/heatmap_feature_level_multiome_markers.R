here::here("metacc/boxplots_feature_level/boxplots_feature_level.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",  help='Cell metadata')
p$add_argument('--met_file',  type="character",  help='DNA methylation file with the feature level quantification')
p$add_argument('--acc_file',  type="character",  help='Chr. accessibility file with the feature level quantification')
p$add_argument('--anno',      type="character",  help='Genomic annotation')
p$add_argument('--markers_file',      type="character",  help='Celltype marker peaks file from the 10x Multiome atlas')
p$add_argument('--outdir',    type="character",  help='Output directory')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

## START TEST ##
args <- list()
args$metadata <- file.path(io$basedir,"results_new/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
args$anno <- "multiome_peaks"
args$met_file <- file.path(io$basedir,sprintf("processed/met/feature_level/%s.tsv.gz",args$anno))
args$acc_file <- file.path(io$basedir,sprintf("processed/acc/feature_level/%s.tsv.gz",args$anno))
args$markers_file <- "/Users/argelagr/data/gastrulation_multiome_10x/results_new/atac/archR/differential/PeakMatrix/markers/marker_peaks_lenient.txt.gz"
args$outdir  <- file.path(io$basedir,"results_new/metacc/boxplots_feature_level/markers")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"absolute"), showWarnings = F)
dir.create(file.path(args$outdir,"relative"), showWarnings = F)

# Options
opts$min_observations <- 50
opts$min_cells <- 10
opts$min_marker_score <- 0.75

# opts$celltypes = c(
#   "Epiblast",
#   "Primitive_Streak",
#   "Caudal_epiblast",
#   "PGC",
#   "Anterior_Primitive_Streak",
#   "Notochord",
#   "Def._endoderm",
#   "Gut",
#   "Nascent_mesoderm",
#   "Mixed_mesoderm",
#   "Intermediate_mesoderm",
#   "Caudal_Mesoderm",
#   "Paraxial_mesoderm",
#   "Somitic_mesoderm",
#   "Pharyngeal_mesoderm",
#   "Cardiomyocytes",
#   "Allantois",
#   "ExE_mesoderm",
#   "Mesenchyme",
#   "Haematoendothelial_progenitors",
#   "Endothelium",
#   "Blood_progenitors",
#   "early_Erythroid",
#   "late_Erythroid",
#   "NMP",
#   "Neurectoderm",
#   "Neural_crest",
#   "Forebrain_Midbrain_Hindbrain",
#   "Spinal_cord",
#   "Surface_ectoderm",
#   "Visceral_endoderm",
#   "ExE_endoderm",
#   "ExE_ectoderm",
#   "Parietal_endoderm"
# )
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

sample_metadata <- fread(args$metadata) %>% 
  .[celltype.mapped%in%opts$celltypes & (pass_metQC==TRUE | pass_accQC==TRUE)]

opts$met_cells <- sample_metadata[pass_metQC==TRUE,id_met]
opts$acc_cells <- sample_metadata[pass_accQC==TRUE,id_acc]

###############
## Load data ##
###############

met.dt <- fread(args$met_file) %>% 
  setnames(c("id_met","id","anno","Nmet","Ntotal","rate")) %>% 
  .[id_met%in%opts$met_cells]

acc.dt <- fread(args$acc_file) %>% 
  setnames(c("id_acc","id","anno","Nmet","Ntotal","rate")) %>% 
  .[id_acc%in%opts$acc_cells]

# Add common cell identifier
met.dt <- merge(met.dt, sample_metadata[,c("cell","id_met")], by="id_met") %>% 
  .[,id_met:=NULL] %>% .[,context:="CG"]
acc.dt <- merge(acc.dt, sample_metadata[,c("cell","id_acc")], by="id_acc") %>% 
  .[,id_acc:=NULL] %>% .[,context:="GC"]

# Merge
metacc.dt <- rbind(met.dt,acc.dt)

rm(met.dt,acc.dt); gc(reset=T)

#######################
## Load marker peaks ##
#######################

markers.to.plot <- c(
  "Epiblast",
  "Primitive_Streak",
  # "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
  "Def._endoderm",
  "Gut",
  # "Nascent_mesoderm",
  "Mixed_mesoderm",
  # "Intermediate_mesoderm",
  # "Caudal_Mesoderm",
  # "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  # "Blood_progenitors",
  "Blood_progenitors_2",
  # "Erythroid1",
  # "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
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
marker_peaks.dt <- fread(args$markers_file) %>%
  .[celltype%in%markers.to.plot]

############################
## Calculate global rates ##
############################

global_rates.dt <- sample_metadata[,c("cell","met_rate","acc_rate")] %>% 
  setnames(c("cell","CG","GC")) %>%
  melt(id.vars="cell", variable.name="context", value.name="global_rate")

##########################
## Plot absolute levels ##
##########################

# celltypes.to.plot <- unique(marker_peaks.dt$celltype)# %>% head(n=3)
# celltypes.to.subset <- sample_metadata[,.N,by=c("celltype.mapped","class")] %>% .[N>=opts$min_cells] %>% .[,.N,by="celltype.mapped"] %>% .[N>1] %>% .$celltype.mapped
  
to.plot <- metacc.dt %>%
  merge(marker_peaks.dt[score>=opts$min_marker_score,c("celltype","idx")] %>% setnames(c("celltype_marker","id")), by="id", allow.cartesian=TRUE) %>%
  merge(sample_metadata[,c("cell","celltype.mapped")], by="cell") %>%
  .[,.(Nmet=sum(Nmet), Ntotal=sum(Ntotal)), by=c("celltype.mapped","celltype_marker","context")] %>% 
  .[,rate:=100*(Nmet/Ntotal)]


ggplot(to.plot[context=="CG"], aes(x = celltype.mapped, y = celltype_marker, fill=rate)) +
  geom_tile(color="black") +
  # facet_wrap(~celltype.mapped) +
  labs(x="", y="") +
  scale_fill_distiller(palette = "YlOrRd", direction=1) +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "top",
    # legend.title = element_blank()
    axis.line = element_blank(),
    axis.text = element_text(color="black"),
    # axis.ticks.x = element_blank()
  )
  

  
pdf(file.path(args$outdir,sprintf("absolute/boxplots_metacc_%s_markers.pdf",i)), width=10, height=6)
print(cowplot::plot_grid(plotlist=list(p.met,p.acc), nrow = 1))
dev.off()


