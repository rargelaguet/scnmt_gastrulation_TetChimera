here::here("metacc/coupling/plot_metacc_coupling.R")

suppressMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Cell metadata')
p$add_argument('--file',    type="character",    help='Precomputed file')
p$add_argument('--outdir',  type="character",    help='Output directory')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
args <- list()
args$metadata <- file.path(io$basedir,"results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
args$file  <- file.path(io$basedir,"results/metacc/coupling/precomputed_metacc_tss_coupling.txt.gz")
args$outdir  <- file.path(io$basedir,"results/metacc/coupling")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"per_cell"), showWarnings = F)
dir.create(file.path(args$outdir,"per_sample"), showWarnings = F)

# Options

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_metQC==TRUE & pass_accQC==TRUE]

# # [TESTING MODE] Subset cells
# if (args$test) {
#   sample_metadata <- sample_metadata %>% .[,head(.SD,n=3),by="sample"]
# }

###########################
## Load precomputed data ##
###########################

metacc_coupling.dt <- fread(args$file)# %>% .[cell%in%sample_metadata$cell]

###################
## Plot per cell ##
###################

cells.to.plot <- unique(metacc_coupling.dt$cell)# %>% head(n=5)

for (i in cells.to.plot) {
  
  to.plot <- metacc_coupling.dt[cell==i]

  p <- ggplot(to.plot, aes(x=window_center, y=r)) +
    stat_summary(fun.data=mean_sd, geom="smooth", alpha=0.2, size=1.0, color="black", fill="black") +
    geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
    geom_vline(xintercept=0, linetype="dashed", color="black", size=0.5) +
    labs(x="Genomic distance from TSS (bp)", y="Met/Acc correlation") +
    coord_cartesian(ylim=c(-0.75,0.5)) +
    theme_classic() +
    theme(
      axis.text = element_text(color="black", size=rel(0.8))
    )
  
  pdf(file.path(args$outdir,sprintf("per_cell/%s.pdf",i)), width=6, height=5)
  print(p)
  dev.off()
}

#####################
## Plot per sample ##
#####################

samples.to.plot <- unique(sample_metadata$sample)

for (i in samples.to.plot) {
  
  to.plot <- metacc_coupling.dt[cell%in%sample_metadata[sample==i,cell]]  %>%
    .[,.(r=mean(r,na.rm=T)), by="window_center"]

  p <- ggplot(to.plot, aes(x=window_center, y=r)) +
    geom_line(size=1.5) +
    # stat_summary(fun.data=mean_sd, geom="smooth", alpha=0.2, size=1.0, color="black", fill="black") +
    geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) +
    geom_vline(xintercept=0, linetype="dashed", color="black", size=0.5) +
    labs(x="Genomic distance from TSS (bp)", y="Met/Acc correlation") +
    coord_cartesian(ylim=c(-0.4,0.05)) +
    theme_classic() +
    theme(
      axis.text = element_text(color="black", size=rel(0.8))
    )

  pdf(file.path(args$outdir,sprintf("per_sample/%s.pdf",i)), width=6, height=5)
  print(p)
  dev.off()
}
