suppressPackageStartupMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--anno',       type="character",  nargs='+',  help='genomic context (i.e. genebody, promoters, etc.')
p$add_argument('--groupA',     type="character",  nargs='+',  help='group A (lineage)')
p$add_argument('--groupB',     type="character",  nargs='+',  help='group B (lineage)')
p$add_argument('--min.cells',  type="integer",                help='Minimum number of cells per group')
p$add_argument('--outfile',    type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))

################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera//settings.R")
  source("/Users/ricard/scnmt_gastrulation_TetChimera//met/differential/utils.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/scnmt_gastrulation_TetChimera//settings.R")
  source("/homes/ricard/scnmt_gastrulation_TetChimera//met/differential/utils.R")
} else {
  stop("Computer not recognised")
}

####################
## Define options ##
####################

# Subset top most variable sites
opts$number_features <- NA

# Filter by coverage
opts$min.CpGs <- 1  # Minimum number of CpG per feature in each cell

# Minimum differential methylation (%) for statistical significance
opts$min.diff <- 25

# Multiple testing correction
opts$threshold_fdr <- 0.10

## START TESTING ##
args$anno <- "ICM_H3K4me3"
args$groupA <- c("hESC")
args$groupB <- c("TE")
args$min.cells <- 10
## END TESTING ##

############################
## Update sample metadata ##
############################

sample_metadata <- sample_metadata %>%
  # .[pass_metQC==TRUE & lineage%in%c(args$groupA,args$groupB)] %>%
  .[!is.na(id_met) & lineage%in%c(args$groupA,args$groupB)] %>%
  .[,group:=as.factor( c("A","B")[as.numeric(lineage%in%args$groupB)+1] )] %>%
  .[,c("id_met","group")]

table(sample_metadata$group)

if (any(sample_metadata[,.N,by="group"]$N<args$min.cells)) {
  stop("Not enough cells per group")
}

###############
## Load data ##
###############

print("Loading data...")


# Load feature metadata
feature_metadata <- fread(
  sprintf("%s/%s.bed.gz",io$features.dir,args$anno), header=FALSE, 
  select=c("V1"="factor", "V2"="integer", "V3"="integer","V5"="character")
) %>% setnames(c("chr","start","end","id"))

# Load methylation data
data <- fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,args$anno), 
  showProgress=F, header=F, select = c("V1"="factor", "V2"= "character", "V4"="integer", "V5"="integer", "V6"="integer")
) %>% .[V1%in%sample_metadata$id_met & V5>=opts$min.CpGs] %>% 
  setnames(c("id_met","id","Nmet","Ntotal","rate"))

# Load methylation data (option 2, cell by cell)
# data <- sample_metadata$id_met %>% map(function(i) {
#   fread(
#     file = sprintf("%s/tmp/%s_%s.gz",io$met_data_parsed,i,args$anno), 
#     showProgress = FALSE, header = FALSE,
#     select = c("V1"="factor", "V2"= "character", "V4"="integer", "V5"="integer", "V6"="integer"))
# }) %>% rbindlist %>% setnames(c("id_met","id","Nmet","Ntotal","rate"))

# Merge methylation data and sample metadata
data <- data %>% merge(sample_metadata, by="id_met")

#############################
## Filter methylation data ##
#############################

print("Filtering data...")

# Filter features by minimum number of cells per group
data <- data[,Ncells:=min(.N), by=c("id","group")] %>% .[Ncells>=args$min.cells] %>% .[,Ncells:=NULL]

# Remove features that have observations in only one group
data <- data[,Ngroup:=length(unique(group)), by="id"] %>% .[Ngroup==2] %>% .[,Ngroup:=NULL]

# Filter by variance
if (!is.na(opts$number_features)) {
  data[,var := var(rate), by="id"] %>% setorder(-var) %>% head(n=opts$number_features) %>% .[,var:=NULL]
}

#########################################
## Differential methylation analysis ##
#########################################

print("Doing differential testing...")

# Binomial assumption: test of equal proportions using Fisher exact test
diff <- data %>% dcast(id~group, value.var=c("Nmet","Ntotal"), fun.aggregate=sum) %>%
  setnames(c("Nmet_A","Nmet_B"),c("A_met","B_met")) %>%
  .[,c("A_unmet","B_unmet"):=list(Ntotal_A-A_met,Ntotal_B-B_met)] %>% 
  .[,c("Ntotal_A","Ntotal_B"):=NULL] %>%
  .[,p.value := fisher.test(x = matrix( c(A_met, A_unmet, B_met, B_unmet), nrow=2, ncol=2))[["p.value"]], by=c("id")] %>%
  .[,c("rateA","rateB"):=list(100*(A_met/(A_met+A_unmet)), 100*(B_met/(B_met+B_unmet)))]

# Multiple testing correction and define significant hits
diff %>%
  .[,diff:=rateB-rateA] %>%
  .[,c("padj_fdr") := list(p.adjust(p.value, method="fdr"))] %>%
  .[,c("log_padj_fdr") := list(-log10(padj_fdr))] %>%
  .[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>opts$min.diff)] %>% 
  .[,(c("rateA","rateB","diff","log_padj_fdr")) := lapply(.SD,round,2), .SDcols = (c("rateA","rateB","diff","log_padj_fdr"))]

# Add genomic coordinates
diff <- merge(feature_metadata, diff, by="id")

# Sort results by p-value
diff %>% setorderv("padj_fdr")

##################
## Save results ##
##################

fwrite(diff, args$outfile, quote=FALSE, sep="\t", na="NA")
