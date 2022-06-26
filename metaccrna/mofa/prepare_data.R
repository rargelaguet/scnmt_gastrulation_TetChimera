###########################
## Load methylation data ##
###########################

met.dt <- lapply(opts$met_annos, function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,n), showProgress=F) %>% 
    .[V1%in%opts$met_cells & V5>=opts$met_min_observations]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","N","rate"))

# Calculate M value from Beta value
met.dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

# Select highly variable features
keep_hv_sites.met <- met.dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id)
met.dt <- met.dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites.met[[y]]]) %>% rbindlist %>% droplevels()

#############################
## Load accessibility data ##
#############################

acc.dt <- lapply(opts$acc_annos, function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$acc_data_parsed,n), showProgress=F) %>% 
    .[V1%in%opts$acc_cells & V5>=opts$acc_min_observations]
}) %>% rbindlist %>% setnames(c("id_acc","id","anno","Nmet","N","rate"))

# Calculate M value from Beta value
acc.dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

# Select highly variable features
keep_hv_sites.acc <- acc.dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% setorder(-var) %>% head(n=opts$acc_nfeatures) %>% .$id)
acc.dt <- acc.dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites.acc[[y]]]) %>% rbindlist %>% droplevels()

###################
## Load RNA data ##
###################

rna.sce <- load_SingleCellExperiment(io$rna.sce, cells=opts$rna_cells, normalise = T)

# Remove lowly variable genes 
# rna.sce <- rna.sce[apply(logcounts(rna.sce),1,var)>=1,]

# Select highly variable features
# keep_hv_genes <- rna.dt %>% .[,.(var=var(expr)), by="ens_id"] %>% setorder(-var) %>% head(n=opts$rna_ngenes) %>% .$ens_id 
# rna.dt <- rna.dt[ens_id%in%keep_hv_genes]
decomp <- modelGeneVar(rna.sce)
hvgs <- decomp[order(decomp$FDR),] %>% head(n=opts$rna_ngenes) %>% rownames

# Convert to data.table
rna.dt <-  as.matrix(logcounts(rna.sce[hvgs,])) %>% t %>% as.data.table(keep.rownames = "id_rna") %>%
  melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id", variable.factor = FALSE)

##############################
## Merge data with metadata ##
##############################

met.dt <- merge(met.dt, cell_metadata.dt[,c("cell","id_met")], by="id_met")
acc.dt <- merge(acc.dt, cell_metadata.dt[,c("cell","id_acc")], by="id_acc")
rna.dt <- merge(rna.dt, cell_metadata.dt[,c("cell","id_rna")], by="id_rna")

###########################
## Prepare data for MOFA ##
###########################

rna.dt <- rna.dt %>% .[,c("cell","ens_id","expr")] %>%
  setnames(c("sample","feature","value")) %>% .[,view:="RNA"]

met.dt <- met.dt %>% .[,c("cell","id","m","anno")] %>%
  setnames(c("sample","feature","value","view")) %>%
  .[,feature:=paste0("met_",feature)] %>%
  .[,view:=paste0("met_",view)]

acc.dt <- acc.dt %>% .[,c("cell","id","m","anno")] %>%
  setnames(c("sample","feature","value","view")) %>%
  .[,feature:=paste0("acc_",feature)] %>%
  .[,view:=paste0("acc_",view)]

data <- do.call("rbind",list(rna.dt,met.dt,acc.dt)) %>%
  .[,value:=round(value,2)]
