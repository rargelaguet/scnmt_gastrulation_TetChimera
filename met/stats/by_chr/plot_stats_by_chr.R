#####################
## Define settings ##
#####################

source("/Users/ricard/scnmt_gastrulation_TetChimera//settings.R")

# Define I/O
io$stats <- paste0(io$basedir,"/met/results/stats/stats_per_chromosome.txt.gz")
io$outdir <- paste0(io$basedir,"/met/results/stats/pdf")

# Define options
# opts$chr <- paste0("chr",c(1:22,"X"))

# Only lineages with enough embryos per sex: 8cell, ICM, TE
opts$lineages <- c(
  # "hESC",
  # "Zygote", 
  # "2cell", 
  # "4cell",
  "8cell",
  # "Morula", 
  "ICM", 
  "TE"
)

# Update sample metadata
sample_metadata <- sample_metadata %>% 
  # .[pass_metQC==TRUE]
  .[!is.na(id_met) & lineage%in%opts$lineages]

sample_metadata[,length(unique(embryo)),by=c("lineage","sex")]

############################
## Load precomputed stats ##
############################

stats <- fread(io$stats) %>%
  .[chr%in%opts$chr] %>% .[,chr:=factor(chr,levels=opts$chr)]

##############################
## Plot global rate per chr ##
##############################

to.plot <- stats %>%
  merge(sample_metadata, by="id_met") %>%
  .[sex%in%c("Female","Male")] %>% .[,sex:=factor(sex,levels=c("Female","Male"))] %>%
  .[,.(rate=mean(mean)), by=c("embryo","chr","sex","lineage")] 

for (i in opts$lineages) {
  p <- ggboxplot(to.plot[lineage==i], x="sex", y="rate", fill="sex") +
    # geom_hline(yintercept=1, linetype="dashed") +
    facet_wrap(~chr, scales="fixed") +
    labs(x="", y="Global DNA methylation rate") +
    theme(
      axis.text.y = element_text(size=rel(0.75)),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  pdf(sprintf("%s/methylation_per_chr_%s.pdf",io$outdir,i), width=8, height=10)
  print(p)
  dev.off()
}
