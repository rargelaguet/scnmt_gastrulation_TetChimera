suppressPackageStartupMessages(library(argparse))

here::here("metacc/qc/qc.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Cell metadata file')
p$add_argument('--outdir',  type="character",              help='Output directory')
p$add_argument('--minimum_number_sites',  type="integer",              help='Minimum number of observations')
p$add_argument('--min_rate',  type="integer",              help='Minimum [met/acc] percentage')
p$add_argument('--max_rate',  type="integer",              help='Maximum [met/acc] percentage')
p$add_argument('--context',  type="character",              help='CG or GC')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$metadata <- file.path(io$basedir,"results/acc/stats/sample_metadata_after_acc_stats.txt.gz")
# args$context <- "GC"
# # args$minimum_number_sites <- 5e3; args$min_rate <- 50; args$max_rate <- 100 # CG
# args$minimum_number_sites <- 1e4; args$min_rate <- 10; args$max_rate <- 40 # GC
# args$outdir <- file.path(io$basedir,"results/acc/qc")
## END TEST ##

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))

# I/O
dir.create(args$outdir, showWarnings = F)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

##############################
## Barplots before QC calls ##
##############################

if (args$context=="CG") {
  to.plot <- sample_metadata %>% .[!is.na(id_met)] %>% setnames(c("nCG","met_rate"),c("N","rate"))
} else if (args$context=="GC") {
  to.plot <- sample_metadata %>% .[!is.na(id_acc)] %>% setnames(c("nGC","acc_rate"),c("N","rate"))
}

to.plot %>% setkey(N) %>% .[,cell:=factor(cell,levels=cell)]

p <- ggplot(to.plot, aes(x=cell, y=log10(N))) +
  geom_bar(stat="identity", position="dodge", fill="gray70", color="gray70") +
  labs(title="", x="", y="Number of observed sites (log10)") +
  geom_hline(yintercept=log10(args$minimum_number_sites), colour="black", linetype="dashed") +
  theme_classic() +
  # scale_y_continuous(expand=c(0,0), limits=c(0,4e+6)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

pdf(file.path(args$outdir,sprintf("%s_qc_barplots_before_qc_calls.pdf",args$context)), width=7, height=5)
print(p)
dev.off()

##############
## QC calls ##
##############

if (args$context=="CG") {
  sample_metadata %>% .[,pass_metQC:=met_rate<=args$max_rate & met_rate>=args$min_rate & nCG>=args$minimum_number_sites]
} else if (args$context=="GC") {
  sample_metadata %>% .[,pass_accQC:=acc_rate<=args$max_rate & acc_rate>=args$min_rate & nGC>=args$minimum_number_sites]
}

#########################################################
## Plot fraction of cells that pass QC for each sample ##
#########################################################

if (args$context=="CG") {
  to.plot <- sample_metadata %>% .[,mean(pass_metQC,na.rm=T), by=c("plate","sample")]
} else if (args$context=="GC") {
  to.plot <- sample_metadata %>% .[,mean(pass_accQC,na.rm=T), by=c("plate","sample")]
}

p <- ggbarplot(to.plot, x="plate", y="V1", fill="gray70") +
  # scale_fill_manual(values=opts$stage.colors) +
  labs(x="", y="Fraction of cells that pass QC") +
  # facet_wrap(~stage)
  theme(
    legend.position = "none",
    axis.text.y = element_text(colour="black",size=rel(0.8)),
    axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
  )

pdf(file.path(args$outdir,sprintf("%s_passQC_barplot.pdf",args$context)), width=7, height=4.5)
print(p)
dev.off()

#############################
## Boxplots after QC calls ##
#############################

if (args$context=="CG") {
  to.plot <- sample_metadata %>% .[pass_metQC==TRUE] %>% setnames(c("nCG","met_rate"),c("N","rate"))
} else if (args$context=="GC") {
  to.plot <- sample_metadata %>% .[pass_accQC==TRUE] %>% setnames(c("nGC","acc_rate"),c("N","rate"))
}

to.plot.melted <- to.plot %>% 
  .[,log10_N:=log10(N)] %>%
  melt(id.vars=c("plate","cell","sample"), measure.vars=c("log10_N","rate")) 

p <- ggplot(to.plot.melted, aes_string(x="plate", y="value")) +
    geom_jitter(size=1, alpha=0.5) +
    geom_boxplot(outlier.shape=NA, coef=1, fill="gray70", alpha=0.8) +
    facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(c("log10_N"="Num. of observations", "rate"="Rate"))) +
    # scale_fill_manual(values=opts$stage.colors) +
    labs(x="", y="") +
    theme_classic() +
    theme(
        axis.text.y = element_text(colour="black",size=rel(1)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
        axis.title.x = element_blank()
    )

pdf(file.path(args$outdir,sprintf("%s_qc_metrics_boxplot_after_qc_calls.pdf",args$context)), width=10, height=5)
print(p)
dev.off()


############################################################
## Scatterplot of mean levels vs number of observed sites ##
############################################################

if (args$context=="CG") {
  to.plot <- copy(sample_metadata) %>% .[!is.na(id_met)] %>% 
    setnames(c("nCG","met_rate"),c("N","rate")) %>% setnames("pass_metQC","pass_QC")
} else if (args$context=="GC") {
  to.plot <- copy(sample_metadata) %>% .[!is.na(id_acc)] %>% 
    setnames(c("nGC","acc_rate"),c("N","rate")) %>% setnames("pass_accQC","pass_QC")
}

p <- ggplot(to.plot, aes(x=rate, y=log10(N), color=pass_QC)) +
  geom_point() +
  geom_hline(yintercept=log10(args$minimum_number_sites), linetype="dashed", color="black") +
  geom_vline(xintercept=args$min_rate, linetype="dashed", color="black") +
  geom_vline(xintercept=args$max_rate, linetype="dashed", color="black") +
  labs(x="Rate (%)", y="Number of observations") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic()

pdf(file.path(args$outdir,sprintf("%s_qc_metrics_scatterplot_rate_vs_N.pdf",args$context)), width=7, height=5)
print(p)
dev.off()

##########
## Save ##
##########

if (args$context=="CG") {
  io$outfile <- file.path(args$outdir, "sample_metadata_after_met_qc.txt.gz")
} else if (args$context=="GC") {
  io$outfile <- file.path(args$outdir, "sample_metadata_after_acc_qc.txt.gz")
}

fwrite(sample_metadata, io$outfile, sep="\t", na = "NA", quote=F)
