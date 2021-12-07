here::here("metacc/profiles/calculate_metacc_profiles.R")

suppressMessages(library(argparse))

######################
## Define arguments ##
######################

# p <- ArgumentParser(description='')
# p$add_argument('--metadata',  type="character",              help='Cell metadata')
# p$add_argument('--indir',  type="character",              help='Input directory')
# p$add_argument('--outdir',  type="character",              help='Output directory')
# p$add_argument('--featuresdir',  type="character",              help='Features directory')
# p$add_argument('--context',  type="character",              help='cg/CG or gc/GC')
# p$add_argument('--annos',    type="character",  nargs="+",  help='Genomic annotation')
# p$add_argument('--test', action="store_true",             help='Test mode? subset number of cells')

# Read arguments
# args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
opts$anno <- "prom_200_200"

# Define window positions and characteristics
opts$positions <- c("prom_200_200"="center")
opts$window_size <- 3000
opts$met.tile <- 100
opts$acc.tile <- 50

opts$test <- TRUE

# I/O
# io$outdir <- file.path(io$basedir,"results/metacc/profiles"); dir.create(io$outdir, showWarnings=F)
io$outfile <- file.path(io$basedir,sprintf("results/metacc/profiles/precomputed_metacc_%s.txt.gz",opts$anno))

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata)

# Define cells
opts$met.cells <- sample_metadata[pass_metQC==TRUE,id_met]
opts$acc.cells <- sample_metadata[pass_accQC==TRUE,id_acc]

if (opts$test) opts$met.cells <- opts$met.cells %>% head(n=5)
if (opts$test) opts$acc.cells <- opts$acc.cells %>% head(n=5)

sample_metadata <- sample_metadata[id_met%in%opts$met_cells | id_acc%in%opts$acc.cells]

##############
## Load data #
##############

# Load genomic annotations
source(here::here("metacc/profiles/load_annotations.R"))

# Load met/acc data
source(here::here("metacc/profiles/load_data.R"))

###########
## Merge ##
###########

# Merge data with sample metadata
met.dt <- met.dt %>% merge(sample_metadata[,c("cell","id_met")], by="id_met") %>% droplevels()
acc.dt <- acc.dt %>% merge(sample_metadata[,c("cell","id_acc")], by="id_acc") %>% droplevels()

# Concatenate DNA methylation and chromatin acessibility data
metacc.dt <- rbind(
  met.dt[,c("cell","id","anno","dist","rate","N","context")],
  acc.dt[,c("cell","id","anno","dist","rate","N","context")]
)

##########
## Save ##
##########

fwrite(metacc.dt, io$outfile)

