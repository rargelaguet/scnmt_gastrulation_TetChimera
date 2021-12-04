suppressPackageStartupMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--anno',       type="character",  nargs='+',  help='genomic context (i.e. genebody, promoters, etc.')
p$add_argument('--lineage',    type="character",  nargs='+',  help='lineage')
p$add_argument('--min.cells',  type="integer",                help='Minimum number of cells per group')
p$add_argument('--outfile',    type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TESTING ##
# args$anno <- "prom_200_200"
# args$lineage <- c("ICM")
# args$min.cells <- 10
## END TESTING ##

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
opts$min.CpGs <- 1                # Minimum number of CpG per feature in each cell

# Minimum differential methylation (%) for statistical significance
opts$min.diff <- 25

# Multiple testing correction
opts$threshold_fdr <- 0.10

############################
## Update sample metadata ##
############################

sample_metadata <- sample_metadata %>%
  .[!is.na(id_met) & lineage%in%args$lineage] %>%
  # .[pass_metQC==TRUE & lineage%in%args$lineage] %>%
  .[,group:=as.factor(sex)] %>%
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

print(object.size(data), units='auto')

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

if (nrow(data)==0) stop("No features passed QC")

#########################################
## Differential methylation analysis ##
#########################################


print("Doing differential testing...")

# Binomial assumption: test of equal proportions using Fisher exact test
diff <- data %>% dcast(id~group, value.var=c("Nmet","Ntotal"), fun.aggregate=sum) %>%
  setnames(c("Nmet_Female","Nmet_Male"),c("Female_met","Male_met")) %>%
  .[,c("Female_unmet","Male_unmet"):=list(Ntotal_Female-Female_met,Ntotal_Male-Male_met)] %>% 
  .[,c("Ntotal_Female","Ntotal_Male"):=NULL] %>%
  .[,p.value := fisher.test(x = matrix( c(Female_met, Female_unmet, Male_met, Male_unmet), nrow=2, ncol=2))[["p.value"]], by=c("id")] %>%
  .[,c("rate_Female","rate_Male"):=list(100*(Female_met/(Female_met+Female_unmet)), 100*(Male_met/(Male_met+Male_unmet)))]
  

# Multiple testing correction and define significant hits
diff %>%
  .[,diff:=rate_Male-rate_Female] %>%
  .[,c("padj_fdr") := list(p.adjust(p.value, method="fdr"))] %>%
  .[,c("log_padj_fdr") := list(-log10(padj_fdr))] %>%
  .[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>opts$min.diff)] %>% 
  .[,(c("rate_Female","rate_Male","diff","log_padj_fdr")) := lapply(.SD,round,2), .SDcols = (c("rate_Female","rate_Male","diff","log_padj_fdr"))]

# Add genomic coordinates
diff <- merge(feature_metadata, diff, by="id")

# Sort results by p-value
diff %>% setorderv("padj_fdr")

##################
## Save results ##
##################

fwrite(diff, args$outfile, quote=FALSE, sep="\t", na="NA")
