#####################
## Define settings ##
#####################

source("/Users/ricard/scnmt_gastrulation_TetChimera/settings.R")

# Define I/O
io$stats <- paste0(io$basedir,"/acc/results/stats/stats_per_chromosome.txt.gz")
io$outdir <- paste0(io$basedir,"/acc/results/stats/pdf")

# Define options
opts$chr <- paste0("chr",c(1:22,"X"))

# Update sample metadata
sample_metadata <- sample_metadata %>% 
  .[pass_accQC==TRUE]
  # .[!is.na(id_acc)]

############################
## Load precomputed stats ##
############################

stats <- fread(io$stats) %>%
  .[chr%in%opts$chr] %>%
  .[,chr:=factor(chr,levels=opts$chr)] %>%
  # .[,mean_coverage:=mean(coverage),by="id_acc"] %>%
  # .[,relative_coverage:=coverage/mean_coverage,by="id_acc"] %>%
  merge(sample_metadata, by="id_acc")

# Sanity check
# all(sort(unique(stats$id_acc)) == sort(sample_metadata$id_acc))

# Normalise by chromosome length
# hg38.genome <- fread(io$hg38.genome) %>%
#   setnames(c("chr","chr_length")) %>%
#   .[chr%in%opts$chr] %>% .[,chr:=factor(chr,levels=opts$chr)]

####################################
## Plot relative coverage per chr ##
####################################

# to.plot <- stats %>%
#   .[,.(coverage=as.double(sum(coverage))), by=c("embryo","chr")] %>%
#   merge(hg38.genome, by="chr") %>% 
#   .[,coverage:=coverage/as.double(chr_length), by="embryo"]# %>%
#   # .[,norm_coverage:=coverage/mean(coverage),by="embryo"]
# 
# ggscatter(to.plot, x="embryo", y="coverage") +
#   # geom_hline(yintercept=1, linetype="dashed") +
#   facet_wrap(~chr, scales="fixed") +
#   labs(x="", y="normalised DNAm coverage") +
#   theme(
#     axis.text.y = element_text(size=rel(0.75)),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )

####################################################################
## Plot accessibility levels per chr, stratifying per sex and lineage ##
####################################################################

to.plot <- stats %>%
  .[sex%in%c("Female","Male")] %>% .[,sex:=factor(sex,levels=c("Female","Male"))] %>%
  .[,.(rate=mean(mean)), by=c("embryo","chr","sex","lineage")] 

sample_metadata[pass_accQC==TRUE,length(unique(embryo)),by=c("lineage","sex")]

for (i in opts$lineages) {
  p <- ggboxplot(to.plot[lineage==i], x="sex", y="rate", fill="sex") +
    # geom_hline(yintercept=1, linetype="dashed") +
    facet_wrap(~chr, scales="fixed") +
    labs(x="", y="Global accessibility rate") +
    theme(
      axis.text.y = element_text(size=rel(0.75)),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  pdf(sprintf("%s/accessibility_per_chr_%s.pdf",io$outdir,i), width=8, height=10)
  # png(sprintf("%s/barplots_%s.png",io$outdir,i), width=6, height=4, units="in", res=400)
  print(p)
  dev.off()
}
