suppressPackageStartupMessages(library(argparse))

here::i_am("metacc/stats/plot_stats.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
# p$add_argument('--query_samples',   type="character",   nargs='+',  help='Query samples')
p$add_argument('--stats',    type="character",  help='Stats file')
p$add_argument('--context',  type="character",              help='cg/CG or gc/GC')
p$add_argument('--celltype_label',  type="character",              help='celltype label')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))

## START TEST ##
args <- list()
args$indir <- file.path(io$basedir,"processed/met/cpg_level")
args$outfile <- file.path(io$basedir,"results/met/stats/sample_met_stats.txt")
args$context <- "CG"
args$celltype_label <- "celltype.mapped_mnn"
## END TEST ##

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

stopifnot(args$celltype_label%in%colnames(sample_metadata))

############################
## Load precomputed stats ##
############################

stats.dt <- fread(args$stats)

# Merge with cell metadata
if (args$context=="CG") {
  stopifnot(stats.dt$id_met%in%sample_metadata$id_met)
  to.plot <- stats %>% merge(sample_metadata, by="id_met") %>% setnames("nCG","N")
} else if (args$context=="GC") {
  stopifnot(stats.dt$id_acc%in%sample_metadata$id_acc)
  to.plot <- stats %>% merge(sample_metadata, by="id_acc") %>% setnames("nGC","N")
}

#######################################
## Boxplots with coverage per sample ##
#######################################

p <- ggboxplot(to.plot, x = "sample", y = "N", outlier.shape=NA) +
  geom_jitter(alpha=0.5, color="#F8766D", size=0.80) +
  yscale("log10", .format = TRUE) +
  labs(x="", y="Number of observations") +
  theme(
    axis.text.y = element_text(size=rel(0.9))
  )

pdf(file.path(args$outdir,sprintf("%s_coverage_per_sample.pdf",args$context)), width=8, height=6)
print(p)
dev.off()

#######################################
## Boxplots with coverage per sample ##
#######################################

p <- ggboxplot(to.plot, x = args$celltype_label, y = "N", outlier.shape=NA) +
  geom_jitter(alpha=0.5, color="#F8766D", size=0.80) +
  yscale("log10", .format = TRUE) +
  labs(x="", y="Number of observations") +
  theme(
    axis.text.y = element_text(size=rel(0.9))
  )

pdf(file.path(args$outdir,sprintf("%s_coverage_per_celltype.pdf",args$context)), width=8, height=6)
print(p)
dev.off()


############################################################
## Correlation between mean methylation rate and coverage ##
############################################################

# to.plot[,mean2:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))]
p <- ggscatter(to.plot, x="rate", y="N") +
  yscale("log10", .format = TRUE) +
  geom_smooth(method="lm") +
  theme_classic()

pdf(file.path(args$outdir,sprintf("%s_correlation_mean_vs_coverage.pdf",args$context)), width=8, height=6)
print(p)
dev.off()
