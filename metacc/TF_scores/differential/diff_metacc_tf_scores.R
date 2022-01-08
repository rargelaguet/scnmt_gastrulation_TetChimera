here::i_am("metacc/differential/feature_level/diff_metacc_feature_level.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Cell metadata')
p$add_argument('--metacc_tf_data',    type="character",  help='Celltypes to use')
p$add_argument('--groupA',     type="character",  nargs='+',  help='group A')
p$add_argument('--groupB',     type="character",  nargs='+',  help='group B')
p$add_argument('--group_label',    type="character",    help='Group label')
p$add_argument('--context',     type="character",  help='CG or GC')
p$add_argument('--min_cells',  type="integer",                help='Minimum number of cells per group')
p$add_argument('--outfile',    type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))


#####################
## Define settings ##
#####################

## START TEST MODE ##
args$metacc_tf_data <- file.path(io$basedir,"results_new/metacc/TF_scores/metacc_tf_binding_sites.tsv.gz")
args$metadata <- file.path(io$basedir,"results_new/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
# args$context <- "GC"
args$groupA <- "KO"
args$groupB <- "WT"
args$group_label <- "ko"
# args$celltype <- "late_Erythroid"
# opts$min_cells <- 10
args$outfile <- file.path(io$basedir,sprintf("results_new/metacc/TF_scores/differential/%s_vs_%s.txt.gz",args$groupA,args$groupB))

## END TEST MODE ##

# Sanity checks
# stopifnot(args$context %in% c("CG","GC"))
# stopifnot(args$celltypes%in%opts$celltypes)

# I/O
dir.create(dirname(args$outfile), showWarnings=F)

####################
## Define options ##
####################

# Define groups
opts$groups <- c(args$groupA,args$groupB)

# Minimum number of observations per TF in each cell
opts$min_features_per_cell <- 25

# Minimum number of cells per celltype and TF
opts$min_cells <- 10

# Statistical test: binomial (counts) or t.test (beta-values)
opts$statistical.test <- "binomial"

# Minimum differential (%) for statistical significance
opts$min.diff <- 25

# Multiple testing correction
opts$threshold_fdr <- 0.10

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>% 
  .[,ko:=ifelse(grepl("KO",class),"KO","WT")] %>%
  .[!is.na(celltype.mapped) & (pass_metQC==TRUE | pass_accQC==TRUE)]

# Define groups
stopifnot(args$group_label%in%colnames(sample_metadata))
sample_metadata <- sample_metadata %>%
  setnames(args$group_label,"group") %>%
  .[group%in%c(args$groupA,args$groupB)] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

# Filter cell types with sufficient number of cells
celltypes.to.use <- sample_metadata[,.N,c("celltype.mapped","group")] %>% .[,sum(N>=opts$min_cells),by="celltype.mapped"] %>% .[V1==2,celltype.mapped]
sample_metadata <- sample_metadata[celltype.mapped%in%celltypes.to.use]
opts$met_cells <- sample_metadata[pass_metQC==TRUE,id_met]
opts$acc_cells <- sample_metadata[pass_accQC==TRUE,id_acc]

# Sanity checks
if (any(sample_metadata[,.N,by=c("group","celltype.mapped")]$N<opts$min_cells)) {
  stop("Some cell types do not have enough cells per group")
}

# print
# cat(args$celltype)
table(sample_metadata$celltype.mapped,sample_metadata$group)

###############
## Load data ##
###############

metacc_tf.dt <- fread(args$metacc_tf_data)

# cat(sprintf("Number of features before filtering:%s \n",length(unique(data$id))))
# cat(sprintf("Number of cells before filtering:%s \n",length(unique(data$cell))))

#################
## Filter data ##
#################

metacc_tf_filt.dt <- metacc_tf.dt %>% .[Ntotal>=opts$min_features_per_cell] %>% 
  merge(sample_metadata[,c("cell","celltype.mapped","group")],by="cell")

# Filter features by minimum number of cells per group
metacc_tf_filt.dt <- metacc_tf_filt.dt %>%
  .[,Ncells:=.N, by=c("tf","group","context")] %>% .[Ncells>=opts$min_cells] %>% .[,Ncells:=NULL]

# cat(sprintf("Number of features after filtering:%s \n",length(unique(metacc_tf_filt.dt$id))))
# cat(sprintf("Number of cells after filtering:%s \n",length(unique(metacc_tf_filt.dt$cell))))

###########################
## Differential analysis ##
###########################

cellypes.to.plot <- unique(metacc_tf_filt.dt$celltype.mapped)

diff.dt <- cellypes.to.plot %>% map(function(i) {
  print(i)
  metacc_tf_filt.dt %>% .[celltype.mapped==i] %>%  
    .[,N:=length(unique(group)), by=c("tf","context")] %>% .[N==2] %>% .[,N:=NULL] %>%  # remove features that have observations in only one group
    dcast(tf+context~group, value.var=c("Nmet","Ntotal"), fun.aggregate=sum) %>%
    setnames(c("tf","context","A_met","B_met","A_total","B_total")) %>%
    .[,c("A_unmet","B_unmet"):=list(A_total-A_met,B_total-B_met)] %>% 
    .[,c("A_total","B_total"):=NULL] %>%
    .[,p.value := fisher.test(x = matrix( c(A_met, A_unmet, B_met, B_unmet), nrow=2, ncol=2))[["p.value"]], by=c("tf","context")] %>%
    .[,c("rateA","rateB"):=list(100*(A_met/(A_met+A_unmet)), 100*(B_met/(B_met+B_unmet)))] %>%
    .[,celltype:=i] %>%
    return
}) %>% rbindlist

# Multiple testing correction and define significant hits
diff.dt %>%
  .[,diff:=rateB-rateA] %>%
  .[,c("padj_fdr") := list(p.adjust(p.value, method="fdr")), by="context"] %>%
  .[,c("log_padj_fdr") := list(-log10(padj_fdr))] %>%
  .[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>opts$min.diff)] %>% 
  .[,(c("rateA","rateB","diff","log_padj_fdr")) := lapply(.SD,round,2), .SDcols = (c("rateA","rateB","diff","log_padj_fdr"))] %>% 
  setorderv("padj_fdr")

##################
## Save results ##
##################

fwrite(diff.dt, args$outfile, quote=FALSE, sep="\t", na="NA")
