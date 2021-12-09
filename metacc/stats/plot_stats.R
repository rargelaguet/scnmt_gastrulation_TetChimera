here::i_am("metacc/stats/plot_stats.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
# p$add_argument('--query_samples',   type="character",   nargs='+',  help='Query samples')
p$add_argument('--metadata',    type="character",  help='Cell metadata file')
p$add_argument('--context',  type="character",              help='cg/CG or gc/GC')
# p$add_argument('--celltype_label',  type="character",              help='celltype label')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$metadata <- file.path(io$basedir,"results/met/stats/sample_metadata_after_met_stats.txt.gz")
# args$outdir <- file.path(io$basedir,"results/met/stats/pdf")
# args$context <- "CG"
# args$celltype_label <- "celltype.mapped_mnn"
## END TEST ##

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))

# I/O
dir.create(args$outdir, showWarnings=F)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

stopifnot(args$celltype_label%in%colnames(sample_metadata))

# Merge with cell metadata
if (args$context=="CG") {
  to.plot <- sample_metadata %>%
    .[!is.na(id_met)] %>%
    .[,c("cell","plate","sample","celltype.mapped_mnn","nCG","met_rate")] %>% 
    setnames(c("nCG","met_rate"),c("N","rate"))
} else if (args$context=="GC") {
  to.plot <- sample_metadata %>%
    .[!is.na(id_acc)] %>%
    .[,c("cell","plate","sample","celltype.mapped_mnn","nGC","acc_rate")] %>% 
    setnames(c("nGC","acc_rate"),c("N","rate"))
}

#######################################
## Boxplots with coverage per sample ##
#######################################

p <- ggboxplot(to.plot, x = "sample", y = "N", outlier.shape=NA, fill="#F8766D", alpha=0.5) +
  geom_jitter(alpha=0.75, fill="#F8766D", size=0.90, shape=21) +
  yscale("log10", .format = TRUE) +
  labs(x="", y="Number of observations") +
  theme(
    axis.text.y = element_text(size=rel(0.80)),
    axis.text.x = element_text(size=rel(0.50))
  )

pdf(file.path(args$outdir,sprintf("%s_coverage_per_sample.pdf",args$context)), width=9, height=4.5)
print(p)
dev.off()

###################################
## Boxplots with rate per sample ##
###################################

p <- ggboxplot(to.plot, x = "sample", y = "rate", outlier.shape=NA, fill="#F8766D", alpha=0.5) +
  geom_jitter(alpha=0.75, fill="#F8766D", size=0.90, shape=21) +
  labs(x="", y="Rate") +
  theme(
    axis.text.y = element_text(size=rel(0.80)),
    axis.text.x = element_text(size=rel(0.50))
  )

pdf(file.path(args$outdir,sprintf("%s_rate_per_sample.pdf",args$context)), width=9, height=4.5)
print(p)
dev.off()


############################################################
## Correlation between mean methylation rate and coverage ##
############################################################

# to.plot[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

p <- ggscatter(to.plot, x="rate", y="N", size=1, color="plate") +
  yscale("log10", .format = TRUE) +
  # xscale("log10", .format = TRUE) +
  geom_smooth(method="lm", color="black") +
  labs(x="Rate", y="Number of observations") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

pdf(file.path(args$outdir,sprintf("%s_correlation_mean_vs_coverage.pdf",args$context)), width=8, height=6)
print(p)
dev.off()
