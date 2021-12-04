#####################
## Define settings ##
#####################

source("/Users/ricard/scnmt_gastrulation_TetChimera/settings.R")

io$outdir <- paste0(io$basedir,"/rna/results/stats")

# update cell metadata
sample_metadata <- sample_metadata[pass_rnaQC==T]

###############
## Load data ##
###############

# SingleCellExperiment object
sce <- readRDS(io$sce)[,sample_metadata$id_rna]

##############################
## Calculate RNA statistics ##
##############################

stats <- perCellQCMetrics(sce) %>% as.data.table %>%
  .[,id_rna:=colnames(sce)] 

##########
## Plot ##
##########

to.plot <- stats %>% merge(sample_metadata, by="id_rna") %>%
  melt(id.vars=c("id_rna","plate"), measure.vars=c("total","detected"), variable.name="metric")

p <- ggboxplot(to.plot, x="plate", y="value", fill="plate") +
  facet_wrap(~metric, scales="free_y") +
  labs(x="", y="") +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

pdf(sprintf("%s/rna_stats.pdf",io$outdir), width=8, height=5, useDingbats = F)
print(p)
dev.off()



to.plot <- stats %>% merge(sample_metadata, by="id_rna") %>%
  .[,.(N=sum(pass_rnaQC)),by=c("plate","tdTOM")]

p <- ggbarplot(to.plot, x="tdTOM", y="N", fill="plate") +
  facet_wrap(~plate) +
  labs(x="", y="Cells that pass RNA QC") +
  theme(
    legend.position = "right",
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

pdf(sprintf("%s/rna_stats.pdf",io$outdir), width=8, height=5, useDingbats = F)
print(p)
dev.off()

##########
## Save ##
##########

fwrite(stats, paste0(io$outdir,"/rna_stats.txt.gz"), sep="\t", quote=F)
