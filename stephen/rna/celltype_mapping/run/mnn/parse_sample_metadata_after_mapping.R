###################
## Load settings ##
###################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else {
  source(here::here("settings.R"))
}

#########
## I/O ##
#########


io$output.metadata <- io$metadata
io$mapping.dir <- paste0(io$basedir,"results/rna/mapping")

#############
## Options ##
#############

# # opts$samples <- c("E8.5_rep1-E8.5_rep2")
# opts$samples <- c(
# 	"E7.5_rep1",
# 	"E7.5_rep2", 
# 	"E8.0_rep1",
# 	"E8.0_rep2",
# 	"E8.5_rep1",
# 	"E8.5_rep2"
# )
opts$cols_to_update <- c("celltype.mapped","celltype.score","closest.cell")

###############
## Load data ##
###############

# Load mapping results
mapping.dt <- dir(io$mapping.dir, full = TRUE, pattern = ".rds") %>%
  map(~readRDS(.x)$mapping %>% setDT) %>%
  rbindlist() %>% 
  .[, .SD, .SDcols = c("cell", opts$cols_to_update)] %>%
  setnames("cell", "id_rna")

# mapping.dt <- opts$samples %>% map(function(x) 
#   readRDS(sprintf("%s/mapping_mnn_%s.rds",io$mapping.dir,x))$mapping %>% .[,c("cell","celltype.mapped","celltype.score","closest.cell")] %>% as.data.table
# ) %>% rbindlist

###########
## Merge ##
###########

sample_metadata <- fread(io$metadata) 
# if columns are already in metadata, remove them first
cols <- colnames(sample_metadata)[colnames(sample_metadata) %in% opts$cols_to_update]
sample_metadata <- sample_metadata[, c(cols) := map(cols, as.null)] %>%
  merge(mapping.dt,by="id_rna",all.x=TRUE)

head(sample_metadata)

#################
## Save output ##
#################

fwrite(sample_metadata, io$output.metadata, sep="\t", na="NA", quote=F)


