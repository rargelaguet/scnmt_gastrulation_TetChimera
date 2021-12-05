here::here("metacc/stats/calculate_stats.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--indir',  type="character",              help='Input directory')
p$add_argument('--outfile',  type="character",              help='Output file')
p$add_argument('--context',  type="character",              help='cg/CG or gc/GC')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
args <- list()
args$indir <- file.path(io$basedir,"processed/acc/gpc_level")
args$outfile <- file.path(io$basedir,"results/acc/stats/sample_acc_stats.txt.gz")
args$context <- "GC"
## END TEST ##

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))

# Define cells
opts$cells <- list.files(args$indir, pattern = "*.tsv.gz") %>% gsub(".tsv.gz","",.)

########################################
## Load data and calculate statistics ##
########################################

stats <- data.table(cell=opts$cells) %>% 
  .[,c("N","rate"):=as.numeric(NA)]

for (i in opts$cells) {
  if (file.exists(sprintf("%s/%s.tsv.gz",args$indir,i))) {
    print(i)

    data <- fread(sprintf("%s/%s.tsv.gz",args$indir,i), sep="\t", verbose=F, showProgress=F, select=c(1,2,4)) %>%
      setnames(c("chr","pos","rate"))

    # Compute genome-wide statistics
    stats[cell==i, c("N","rate"):=list(nrow(data), round(100*mean(data$rate),2))]

  } else {
    print(sprintf("Sample %s not found",i))
  }
}

# Define column names
if (args$context=="CG") {
  stats %>% setnames(c("cell","nCG","rate"))
} else if (args$context=="GC") {
  stats %>% setnames(c("cell","nGC","rate"))
}

##########
## Save ##
##########

fwrite(stats, args$outfile, sep="\t", na = "NA", quote=F)
