suppressMessages(library(MOFA2))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera/acc/dimensionality_reduction/load_settings.R")
} else {
  stop("Computer not recognised")
}

#############################
## Load accessibility data ##
#############################

acc_dt <- lapply(names(opts$annos), function(n)
  fread(sprintf("%s/%s.tsv.gz",io$acc_data_parsed,n), select=c(1,2,3,5,6)) %>% 
    setnames(c("id_acc","id","anno","N","rate")) %>% 
    .[id_acc%in%sample_metadata$id_acc] %>%
    .[N>=opts$min.GpCs] %>% .[,N:=NULL]
) %>% rbindlist

# Merge with sample metadata 
acc_dt <- acc_dt %>% merge(sample_metadata, by="id_acc")

#######################################
## Calculate M value from Beta value ##
#######################################

acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

#################
## Filter data ##
#################

# Filter features by coverage
nsamples <- length(unique(acc_dt$id_acc))
acc_dt <- acc_dt[,cov:=.N/nsamples,by=c("id","anno")] %>% .[cov>=opts$min.coverage] %>% .[,c("cov"):=NULL]

# Regress out global accessibility rate
# foo <- acc_dt[,.(mean=mean(m)),by=c("id_acc","anno")]
# acc_dt <- acc_dt %>% merge(foo, by=c("id_acc","anno")) %>%
#   .[,m:=mean(m)+lm(formula=m~mean)[["residuals"]], by=c("id","anno","stage_lineage")]

# Filter features by variance
keep_hv_sites <- acc_dt %>% split(.$anno) %>% map(~ .[,.(var = var(rate)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n = opts$nfeatures) %>% .$id)
acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist

###########################
## Prepare data for MOFA ##
###########################

data <- acc_dt %>% .[,c("id_acc","id","m","anno")] %>%
  setnames(c("sample","feature","value","view"))

##########
## Save ##
##########

# file <- paste0(io$outdir,"/data.txt")
# fwrite(data, file, col.names=T, quote=F, sep="\t")
# system(sprintf("pigz -f %s", file))


####################
## Fit MOFA model ##
####################

MOFAobject <- create_mofa(data)

# Visualise data structure
plot_data_overview(MOFAobject)

####################
## Define options ##
####################

# Data options
data_opts <- get_default_data_options(MOFAobject)

# model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 3

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42

#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# Train the model
model <- run_mofa(MOFAobject, io$outfile)
