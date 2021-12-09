here::i_am("metacc/differential/feature_level/diffmet_feature_level.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
# p$add_argument('--metadata',  type="character",              help='Cell metadata')
p$add_argument('--anno',       type="character",  help='genomic context (i.e. genebody, promoters, etc.')
p$add_argument('--groupA',     type="character",  nargs='+',  help='group A')
p$add_argument('--groupB',     type="character",  nargs='+',  help='group B')
p$add_argument('--group_label',    type="character",    help='Group label')
p$add_argument('--context',     type="character",  nargs='+',  help='CG or GC')
p$add_argument('--min_cells',  type="integer",                help='Minimum number of cells per group')
p$add_argument('--outfile',    type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TESTING ##
# args$metadata <- file.path(io$basedir,"results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
# args$anno <- "prom_200_200"
# args$groupA <- c("Surface_ectoderm_E8.5_chimera_KO")
# args$groupB <- c("Surface_ectoderm_E8.5_host_WT")
# args$group_label <- "celltype_class"
# args$context <- "CG"
# args$min_cells <- 10
# args$outfile <- file.path(io$basedir,sprintf("results/met/differential/%s/%s_%s_vs_%s.txt.gz",args$group_label,args$group_label,args$groupA,args$groupB))
## END TESTING ##

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))

# I/O
dir.create(dirname(args$outfile))

####################
## Define options ##
####################

# Define groups
opts$groups <- c(args$groupA,args$groupB)

# Minimum number of CpG per feature in each cell
opts$min_obs <- 1                

# Statistical test: binomial (counts) or t.test (beta-values)
opts$statistical.test <- "binomial"

# Minimum differential (%) for statistical significance
opts$min.diff <- 25

# Multiple testing correction
opts$threshold_fdr <- 0.10

##########################
## Load sample metadata ##
##########################

print("Loading cell metadata...")

sample_metadata <- fread(io$metadata) %>%
  .[,celltype_class:=sprintf("%s_%s",celltype.mapped_mnn,class)] %>% # temporary
  .[pass_metQC==TRUE]

stopifnot(args$group_label%in%colnames(sample_metadata))

sample_metadata <- sample_metadata %>%
  setnames(args$group_label,"group") %>%
  # .[,group:=eval(as.name(args$group_label))] %>%
  .[group%in%c(args$groupA,args$groupB)] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

# test mode
if (isTRUE(args$test_mode)) {
  print("Testing mode activated")
  sample_metadata <- sample_metadata %>% split(.,.$group) %>% map(~ head(.,n=250)) %>% rbindlist
} else {
  # Sanity checks
  if (any(sample_metadata[,.N,by="group"]$N<args$min_cells)) {
    stop("Not enough cells per group")
  }
}

table(sample_metadata$group)

# Define cells to use
if (args$context=="CG") {
  cells <- sample_metadata$id_met
  io$indir <- io$met_data_parsed
  sample_metadata <- sample_metadata[,c("id_met","cell","group")]
} else {
  cells <- sample_metadata$id_acc
  io$indir <- io$acc_data_parsed
  sample_metadata <- sample_metadata[,c("id_acc","cell","group")]
}
###########################
## Load feature metadata ##
###########################

print("Loading feature metadata...")

feature_metadata <- fread(
  file.path(io$features.dir,sprintf("%s.bed.gz",args$anno)), header=FALSE, 
  select=c("V1"="factor", "V2"="integer", "V3"="integer","V5"="character")
) %>% setnames(c("chr","start","end","id"))

###############
## Load data ##
###############

print("Loading data...")

# Load methylation data
data <- fread(file.path(io$indir,sprintf("%s.tsv.gz",args$anno)), 
  showProgress=F, header=F, select = c("V1"="factor", "V2"= "character", "V4"="integer", "V5"="integer", "V6"="integer")
) %>% .[V1%in%cells & V5>=opts$min_obs] %>% 
  setnames(c("id_met","id","Nmet","Ntotal","rate"))

# Load methylation data (option 2, cell by cell)
# data <- sample_metadata$id_met %>% map(function(i) {
#   fread(
#     file = sprintf("%s/tmp/%s_%s.gz",io$met_data_parsed,i,args$anno), 
#     showProgress = FALSE, header = FALSE,
#     select = c("V1"="factor", "V2"= "character", "V4"="integer", "V5"="integer", "V6"="integer"))
# }) %>% rbindlist %>% setnames(c("id_met","id","Nmet","Ntotal","rate"))

print(object.size(data), units='auto')

# Merge methylation data and sample metadata
data <- data %>% merge(sample_metadata, by="id_met")

# Convert beta value to M value
if (opts$statistical.test == "t.test") {
  data[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]
}

#############################
## Filter methylation data ##
#############################

print("Filtering data...")

# Filter features by minimum number of cells per group
data <- data[,Ncells:=min(.N), by=c("id","group")] %>% .[Ncells>=args$min_cells] %>% .[,Ncells:=NULL]

# Remove features that have observations in only one group
data <- data[,Ngroup:=length(unique(group)), by="id"] %>% .[Ngroup==2] %>% .[,Ngroup:=NULL]

# Filter by variance
# if (!is.na(opts$number_features)) {
#   data[,var := var(rate), by="id"] %>% setorder(-var) %>% head(n=opts$number_features) %>% .[,var:=NULL]
# }

#########################################
## Differential methylation analysis ##
#########################################

print("Doing differential testing...")

# Binomial assumption: test of equal proportions using Fisher exact test
if (opts$statistical.test == "binomial") {
  diff.dt <- data %>% 
    dcast(id~group, value.var=c("Nmet","Ntotal"), fun.aggregate=sum) %>%
    setnames(c("id","A_met","B_met","A_total","B_total")) %>%
    .[,c("A_unmet","B_unmet"):=list(A_total-A_met,B_total-B_met)] %>% 
    .[,c("A_total","B_total"):=NULL] %>%
    .[,p.value := fisher.test(x = matrix( c(A_met, A_unmet, B_met, B_unmet), nrow=2, ncol=2))[["p.value"]], by=c("id")] %>%
    .[,c("rateA","rateB"):=list(100*(A_met/(A_met+A_unmet)), 100*(B_met/(B_met+B_unmet)))]
  
# T-test under normality assumption
} else if (opts$statistical.test == "t.test") {
  stop("Not implemented")
}

# Multiple testing correction and define significant hits
diff.dt %>%
  .[,diff:=rateB-rateA] %>%
  .[,c("padj_fdr") := list(p.adjust(p.value, method="fdr"))] %>%
  .[,c("log_padj_fdr") := list(-log10(padj_fdr))] %>%
  .[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>opts$min.diff)] %>% 
  .[,(c("rateA","rateB","diff","log_padj_fdr")) := lapply(.SD,round,2), .SDcols = (c("rateA","rateB","diff","log_padj_fdr"))]

# Add genomic coordinates
# diff.dt <- merge(feature_metadata, diff.dt, by="id")

# Sort results by p-value
diff.dt %>% setorderv("padj_fdr")

##################
## Save results ##
##################

fwrite(diff.dt, args$outfile, quote=FALSE, sep="\t", na="NA")
