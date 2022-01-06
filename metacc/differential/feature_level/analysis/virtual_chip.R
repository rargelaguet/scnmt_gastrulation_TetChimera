here::i_am("metacc/differential/feature_level/analysis/virtual_chip.R")

#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
opts$motif_annotation <- "CISBP"

# opts$TFsTFs <- colnames(virtual_chip.mtx)
opts$TFs <- c("TAL1", "GATA1", "RUNX1", "FOXA2", "GATA4", "CDX2","NKX2-5","TBX5", "SOX10")

# I/O
io$virtual_chip.dir <- file.path(io$multiome.basedir,"results_new/rna_atac/virtual_chipseq/CISBP")
io$virtual_chip.mtx <- file.path(io$virtual_chip.dir,"virtual_chip.mtx")
io$outdir <- file.path(io$basedir,"results_new/met/differential/test"); dir.create(io$outdir, showWarnings = F)

##################################################
## Load differential results for multiome peaks ##
##################################################

diff.dt <- fread("/Users/argelagr/data/tet_chimera_nmtseq/results_new/met/differential/multiome_peaks_Blood_progenitors_WT_vs_KO.txt.gz")

###################################
## Load virtual ChIP-seq library ##
###################################

# Load detailed data.tables
# virtual_chip.dt <- list.files(io$virtual_chip.dir,".txt.gz") %>% strsplit("\\.") %>% map_chr(1) %>% map(function(i) {
#     fread(sprintf("%s/%s.txt.gz",io$virtual_chip.dir,i)) %>%
#     # .[,c("chr","start","end"):=NULL] %>%
#     .[,tf:=i] %>%
#     return
# }) %>% rbindlist

# Load matrix
virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)

#################
## Filter data ##
#################

peaks <- intersect(diff.dt$id,rownames(virtual_chip.mtx))
diff.dt <- diff.dt[id%in%peaks]
virtual_chip.mtx <- virtual_chip.mtx[peaks,]

virtual_chip.mtx <- virtual_chip.mtx[,colSums(virtual_chip.mtx)>=25]

#########
## foo ##
#########

# plot differential methylation for each group of peaks that are bound by the same TF
i <- "GATA1"
diff_tf.dt <- colnames(virtual_chip.mtx) %>% map(function(i) {
  tf.peaks <- names(which(virtual_chip.mtx[,i]>=0.25))
  tmp <- diff.dt[id%in%tf.peaks] %>% .[,tf:=i]
  
  nrow(tmp[sig==T])
  length(tf.peaks)
  
  nrow(tmp[sig==F])
  nrow(tmp)
  
  fisher.test()
  fisher.test(x = matrix( c(A_met, A_unmet, B_met, B_unmet), nrow=2, ncol=2))[["p.value"]]
    
  return(tmp)
}) %>% rbindlist

diff_tf.dt[,.(fraction=mean(sig==T), N=sum(sig==T)),by=c("tf")] %>% .[N>=25] %>% View

diff_tf.dt[sig==T & diff>0,.N,by="tf"] %>% View
diff_tf.dt[tf=="ETV2"]

tf.peaks
# Hypergeometric test

###########
## Plot  ##
###########

ggboxplot()


p <- ggplot(to.plot[N>0], aes_string(x="min_score", y="log10_N", color="cor_sign")) +
  geom_line(size=1) +
  facet_wrap(~tf, scales="free_y") +
  labs(y="Number of predicted binding sites (log10)", x="Minimum in silico binding score") +
  scale_color_brewer(palette="Dark2") +
  # coord_cartesian(xlim=c(0,0.9)) +
  theme_classic() +
  theme(
    axis.text = element_text(color="black"),
    # legend.position = "right",
    legend.position = "top",
    legend.title = element_blank()
  )

pdf(file.path(io$outdir,"lineplot_number_binding_sites_per_tf_positive_vs_negative.pdf"), width=6, height=5)
print(p)
dev.off()

