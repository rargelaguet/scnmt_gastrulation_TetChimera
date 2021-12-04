suppressPackageStartupMessages(library(argparse))

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('--anno',       type="character",              help='genomic context (i.e. genebody, promoters, etc.')
p$add_argument('--groupA',     type="character",  nargs='+',  help='group 1 (lineage)')
p$add_argument('--groupB',     type="character",  nargs='+',  help='group 2 (lineage)')
p$add_argument('--outfile',    type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TESTING ##
# args$anno <- "LTR"
# args$groupA <- c("ICM")
# args$groupB <- c("TE_polar")
## END TESTING ##

###########################
## Load default settings ##
###########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera//settings.R")
  # source("/Users/ricard/scnmt_gastrulation_TetChimera/met/differential/utils.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/scnmt_gastrulation_TetChimera//settings.R")
  # source("/homes/ricard/scnmt_gastrulation_TetChimera/met/differential/utils.R")
} else {
  stop("Computer not recognised")
}
io$outfile <- args$outfile
io$met_data_parsed <- paste0(io$met_data_parsed,"/pseudobulk")

opts$min.observations <- 5

###############
## Load data ##
###############

# Load methylation data
data <- fread(
  file = sprintf("%s/%s.tsv.gz",io$met_data_parsed,args$anno), 
  showProgress=F, header=F,
  select = c("V1"="factor", "V2"= "character", "V4"="integer", "V5"="integer", "V6"="integer")
) %>% .[V1%in%c(args$groupA,args$groupB)] %>% droplevels %>% setnames(c("lineage","id","Nmet","N","rate"))

# Load genomic coordinates
feature_metadata <- fread(
  sprintf("%s/%s.bed.gz",io$features.dir,args$anno), header=FALSE, 
  select=c("V1"="factor", "V2"="integer", "V3"="integer","V5"="character")
) %>% setnames(c("chr","start","end","id"))

############################
## Parse methylation data ##
############################

diff <- data %>% 
  .[N>=opts$min.observations] %>%
  .[lineage==args$groupA,lineage:="A"] %>% .[lineage==args$groupB,lineage:="B"] %>%
  dcast(id~lineage, value.var=c("Nmet","N","rate")) %>%
  .[,diff_rate:=rate_B-rate_A] %>%
  .[!is.na(diff_rate)] 
  
# Add genomic coordinates
diff <- merge(feature_metadata,diff, by="id")

##################
## Save results ##
##################

fwrite(diff, io$outfile, quote=FALSE, sep="\t", na="NA")
