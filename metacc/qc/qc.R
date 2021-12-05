suppressPackageStartupMessages(library(argparse))

here::here("metacc/stats/calculate_stats.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
# p$add_argument('--indir',  type="character",              help='Input directory')
# p$add_argument('--stats',    type="character",  help='Stats file')
p$add_argument('--metadata',  type="character",              help='Cell metadata file')
p$add_argument('--outfile',  type="character",              help='Output file')
p$add_argument('--minimum_number_sites',  type="integer",              help='Minimum number of observations')
p$add_argument('--min_rate',  type="integer",              help='Minimum [met/acc] percentage')
p$add_argument('--max_rate',  type="integer",              help='Maximum [met/acc] percentage')
p$add_argument('--context',  type="character",              help='CG or GC')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
args <- list()
args$stats <- file.path(io$basedir,"results/met/stats/sample_met_stats.txt.gz")
args$metadata <- file.path(io$basedir,"results/met/stats/sample_metadata_after_met_stats.txt.gz")
args$context <- "CG"
args$minimum_number_sites <- 1e4
args$min_rate <- 10
args$max_rate <- 50
args$outfile <- file.path(io$basedir,"results/met/qc/sample_metadata_after_met_qc.txt.gz")
## END TEST ##

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))

# Define cells
# opts$cells <- list.files(args$indir, pattern = "*.tsv.gz") %>% gsub(".tsv.gz","",.)

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

##############################
## Barplots before QC calls ##
##############################

tmp <- to.plot[,c("id_met","nCG")] %>% setkey(nCG) %>% .[,id_met:=factor(id_met,levels=id_met)]

p <- ggplot(tmp, aes(x=id_met, y=log10(nCG))) +
  geom_bar(stat="identity", position="dodge", fill="#F8766D", color="#F8766D") +
  labs(title="", x="", y="Number of observed CpG sites (log10)") +
  geom_hline(yintercept=log10(opts$met_nCG_threshold), colour="black", linetype="dashed") +
  barplot_theme() +
  # scale_y_continuous(expand=c(0,0), limits=c(0,4e+6)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
print(p)

pdf(file=file.path(args$outdir,"qc_met.pdf"), width=8, height=5)
print(p1)
dev.off()

##############
## QC calls ##
##############

sample_metadata %>% .[,pass_metQC:=met_rate<=args$max_rate & met_rate>=args$min_rate & met_sites>=args$minimum_number_sites]


#############################
## Boxplots after QC calls ##
#############################

to.plot <- metadata %>% .[pass_rnaQC==TRUE] %>%
    melt(id.vars=c("plate","cell","stage"), measure.vars=c("nFeature_RNA","mit_percent_RNA","rib_percent_RNA"))

facet.labels <- c("nFeature_RNA" = "Num. of genes", "mit_percent_RNA" = "Mitochondrial %", "rib_percent_RNA" = "Ribosomal %")
    

p <- ggplot(to.plot, aes_string(x="plate", y="value", fill="stage")) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    geom_jitter(shape=21) +
    facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(facet.labels)) +
    # scale_fill_manual(values=opts$stage.colors) +
    theme_classic() +
    theme(
        axis.text.y = element_text(colour="black",size=rel(1)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
        axis.title.x = element_blank()
    )

pdf(file.path(args$outdir,"qc_metrics_boxplot.pdf"), width=9, height=5)
# pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outdir))
print(p)
dev.off()




print("Fail QC for methylation:")
failqc <- stats %>%
  # .[((rate<opts$min.meth.rate) | (rate>opts$max.meth.rate))] %>%
  .[(nCG<opts$met_nCG_threshold) | (rate<opts$min.meth.rate | rate>opts$max.meth.rate)] %>%
  .[,id_met]

print(length(failqc))


############################################################
## Scatterplot of mean levels vs number of observed sites ##
############################################################

to.plot <- stats %>%
  .[,pass_metQC:=ifelse(id_met%in%failqc,F,T)]

p <- ggplot(to.plot, aes(x=rate, y=log10(nCG), color=pass_metQC)) +
  geom_point() +
  geom_hline(yintercept=log10(opts$met_nCG_threshold), linetype="dashed", color="black") +
  geom_vline(xintercept=opts$min.meth.rate, linetype="dashed", color="black") +
  geom_vline(xintercept=opts$max.meth.rate, linetype="dashed", color="black") +
  labs(x="Global methylation (%)", y="Number of observed CpG sites") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic()



##########
## Save ##
##########



fwrite(stats, args$outfile, sep="\t", na = "NA", quote=F)
