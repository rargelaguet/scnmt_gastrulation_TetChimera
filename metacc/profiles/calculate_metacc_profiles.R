here::here("metacc/profiles/calculate_metacc_profiles.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Cell metadata')
p$add_argument('--anno',    type="character",  help='Genomic annotation')
p$add_argument('--window_size',  type="integer", default=3000,              help='Window size')
p$add_argument('--met_tile',  type="integer", default=100,              help='DNA methylation tile size')
p$add_argument('--acc_tile',  type="integer", default=100,              help='Chromatin accessibility tile size')
p$add_argument('--test', action="store_true",             help='Test mode? subset number of cells')
p$add_argument('--outfile',  type="character",              help='Output file')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Define settings ##
###################

## START TEST ##
# args <- list()
# args$metadata <- file.path(io$basedir,"results_new/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
# args$anno <- "multiome_peaks"
# args$window_size <- 1500
# args$met_tile <- 150
# args$acc_tile <- 75
# args$outfile  <- file.path(io$basedir,sprintf("results_new/metacc/profiles/%s/precomputed_metacc_%s_v2.txt.gz",args$anno,args$anno))
# args$test <- TRUE
## END TEST ##

# I/O
dir.create(dirname(args$outfile), showWarnings = F)

# Options

# Define window positions and characteristics
opts$positions <- c("center"); names(opts$positions) <- args$anno

###################
## Load metadata ##
###################

sample_metadata <- fread(args$metadata)

# Define cells
opts$met.cells <- sample_metadata[pass_metQC==TRUE,id_met]
opts$acc.cells <- sample_metadata[pass_accQC==TRUE,id_acc]
# opts$met.cells <- sample_metadata[!is.na(id_met),id_met]
# opts$acc.cells <- sample_metadata[!is.na(id_acc),id_acc]

if (args$test) opts$met.cells <- opts$met.cells %>% head(n=5)
if (args$test) opts$acc.cells <- opts$acc.cells %>% head(n=5)

sample_metadata <- sample_metadata[id_met%in%opts$met.cells | id_acc%in%opts$acc.cells]

table(sample_metadata$sample)

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
print(sprintf("DNA methylation, number of cells before merging: %s",length(unique(met.dt$id_met))))
met.dt <- met.dt %>% merge(sample_metadata[,c("cell","id_met")], by="id_met") %>% droplevels()
print(sprintf("DNA methylation, number of cells after merging: %s",length(unique(met.dt$cell))))
print(sprintf("Cells not present in the DNA methylation data: %s",paste(opts$met.cells[!opts$met.cells%in%unique(met.dt$id_met)], collapse=" ")))

print(sprintf("Chromatin accessibility, number of cells before merging: %s",length(unique(acc.dt$id_acc))))
acc.dt <- acc.dt %>% merge(sample_metadata[,c("cell","id_acc")], by="id_acc") %>% droplevels()
print(sprintf("Chromatin accessibility, number of cells after merging: %s",length(unique(acc.dt$cell))))
print(sprintf("Cells not present in the chromatin accessibility data: %s",paste(opts$acc.cells[!opts$acc.cells%in%unique(acc.dt$id_acc)], collapse=" ")))

# Concatenate DNA methylation and chromatin acessibility data
metacc.dt <- rbind(
  met.dt[,c("cell","id","anno","dist","rate","N","context")],
  acc.dt[,c("cell","id","anno","dist","rate","N","context")]
)


##########
## Save ##
##########

fwrite(metacc.dt, args$outfile)

