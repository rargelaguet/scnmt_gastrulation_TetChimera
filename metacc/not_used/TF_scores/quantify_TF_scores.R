# here::i_am("")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$met_file <- file.path(io$basedir,"processed/met/feature_level/multiome_peaks.tsv.gz")
io$acc_file <- file.path(io$basedir,"processed/acc/feature_level/multiome_peaks.tsv.gz")
io$virtual_chip.dir <- file.path(io$multiome.basedir,"results_new/rna_atac/virtual_chipseq/CISBP")
io$virtual_chip.mtx <- file.path(io$virtual_chip.dir,"virtual_chip.mtx")
io$outdir <- file.path(io$basedir,"results_new/metacc/TF_scores"); dir.create(io$outdir, showWarnings = F)

# Options
opts$motif_annotation <- "CISBP"

# opts$TFsTFs <- colnames(virtual_chip.mtx)
opts$TFs <- c("TAL1", "GATA1", "RUNX1", "FOXA2", "GATA4", "CDX2","NKX2-5","TBX5", "SOX10")

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>% 
  .[!is.na(celltype.mapped) & (pass_metQC==TRUE | pass_accQC==TRUE)]

opts$met_cells <- sample_metadata[pass_metQC==TRUE,id_met]
opts$acc_cells <- sample_metadata[pass_accQC==TRUE,id_acc]

####################################################
## Load met/acc quantification for multiome peaks ##
####################################################

met.dt <- fread(io$met_file) %>% 
  setnames(c("id_met","id","anno","Nmet","Ntotal","rate")) %>% 
  .[id_met%in%opts$met_cells]

acc.dt <- fread(io$acc_file) %>% 
  setnames(c("id_acc","id","anno","Nmet","Ntotal","rate")) %>% 
  .[id_acc%in%opts$acc_cells]

# Add common cell identifier
met.dt <- merge(met.dt, sample_metadata[,c("cell","id_met")], by="id_met") %>% 
  .[,id_met:=NULL] %>% .[,context:="CG"]
acc.dt <- merge(acc.dt, sample_metadata[,c("cell","id_acc")], by="id_acc") %>% 
  .[,id_acc:=NULL] %>% .[,context:="GC"]

# Merge
metacc.dt <- rbind(met.dt,acc.dt)

###################################
## Load virtual ChIP-seq library ##
###################################

virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)

#################
## Filter data ##
#################

peaks <- intersect(unique(metacc.dt$id),rownames(virtual_chip.mtx))
virtual_chip.mtx <- virtual_chip.mtx[peaks,]
metacc.dt <- metacc.dt[id%in%peaks,]

############################################
## Quantify met/acc over TF binding sites ##
############################################

opts$min.chip.score <- 0.25

TFs.to.plot <- colnames(virtual_chip.mtx)

# i <- "GATA1"
metacc_tf.dt <- TFs.to.plot %>% map(function(i) {
  tf.peaks <- names(which(virtual_chip.mtx[,i]>=opts$min.chip.score))
  metacc.dt[id%in%tf.peaks] %>% 
    .[,.(Nmet=sum(Nmet), Ntotal=sum(Ntotal)), by=c("cell","context")] %>% 
    .[,rate:=round(100*(Nmet/Ntotal))] %>%
    .[,tf:=i] %>%
    return
}) %>% rbindlist

# Save
fwrite(metacc_tf.dt, file.path(io$outdir,"metacc_tf_binding_sites.tsv.gz"), sep="\t")

############################
## Calculate global rates ##
############################

global_rates.dt <- metacc.dt %>%
  .[,.(Nmet=sum(Nmet), Ntotal=sum(Ntotal)), by=c("cell","context")] %>% 
  .[,global_rate:=round(100*(Nmet/Ntotal))] %>% .[,c("Nmet","Ntotal"):=NULL]
  
# global_rates.dt <- sample_metadata[,c("cell","met_rate","acc_rate")] %>% 
#   setnames(c("cell","CG","GC")) %>%
#   melt(id.vars="cell", variable.name="context", value.name="global_rate")

##########
## Plot ##
##########

to.plot <- metacc_tf.dt %>% .[Ntotal>=15] %>% 
  merge(sample_metadata,by="cell") %>%
  merge(global_rates.dt, by=c("cell","context")) %>% 
  .[,rate_norm:=log2(rate/global_rate)]

TFs.to.plot <- unique(to.plot$tf)

for (i in TFs.to.plot) {
  
  celltypes.to.plot.met <- to.plot[context=="CG" & tf==i] %>% 
    .[,.(N=.N),by=c("celltype.mapped","class")] %>% .[,sum(N>=5),by="celltype.mapped"] %>% .[V1==2,celltype.mapped]
  celltypes.to.plot.acc <- to.plot[context=="GC" & tf==i] %>% 
    .[,.(N=.N),by=c("celltype.mapped","class")] %>% .[,sum(N>=5),by="celltype.mapped"] %>% .[V1==2,celltype.mapped]
  celltypes.to.plot <- intersect(celltypes.to.plot.met,celltypes.to.plot.acc)
  
  if (length(celltypes.to.plot)>=1) {
    
    to.plot2 <- to.plot[context=="CG" & tf==i & celltype.mapped%in%celltypes.to.plot]
    
    p.met <- ggplot(to.plot2, aes(x = class, y = rate_norm, fill=class)) +
      geom_boxplot(outlier.shape=NA, coef=1) +
      geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 1, alpha = 0.5, shape=21) +
      facet_wrap(~celltype.mapped) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(x="", y=sprintf("%s Methylation levels (%%)",i), title=i) +
      theme_bw() +
      guides(x = guide_axis(angle = 90)) +
      theme(
        plot.title = element_text(hjust = 0.5, size=rel(1)),
        legend.position = "top",
        legend.title = element_blank(),
        # axis.text = element_text(color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    
    to.plot2 <- to.plot[context=="GC" & tf==i & celltype.mapped%in%celltypes.to.plot]
    
    p.acc <- ggplot(to.plot2, aes(x = class, y = rate_norm, fill=class)) +
      geom_boxplot(outlier.shape=NA, coef=1) +
      geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 1, alpha = 0.5, shape=21) +
      facet_wrap(~celltype.mapped) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(x="", y=sprintf("%s Accessibility levels (%%)",i), title=i) +
      theme_bw() +
      guides(x = guide_axis(angle = 90)) +
      theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size=rel(1)),
        legend.title = element_blank(),
        # axis.text = element_text(color="black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    
    pdf(file.path(io$outdir,sprintf("%s_boxplots_metacc.pdf",i)), width=10, height=6)
    print(cowplot::plot_grid(plotlist=list(p.met,p.acc), nrow = 1))
    dev.off()
  }
}
