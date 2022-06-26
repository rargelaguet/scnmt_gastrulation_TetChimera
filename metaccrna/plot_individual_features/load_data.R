scale <- function(value, min.data, max.data, min.scaled, max.scaled) {
  stopifnot(is.numeric(value))
  stopifnot(value<=max.data & value>=min.data)
  return ((max.scaled - min.scaled) * (value - min.data) / (max.data - min.data)) + min.scaled
}

load_data <- function(io, rna_gene, met_feature, met_context, acc_feature, acc_context, min_cpg=1, min_gpc=1) {
  
  # Load DNA methylation data
  met_dt <- fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,met_context)) %>%
    setnames(c("id_met","id","anno","Nmet","N","value")) %>%
    .[id%in%met_feature] %>% .[N>=min_cpg]
  
  # Load DNA accessibility data
  acc_dt <- fread(sprintf("%s/%s.tsv.gz",io$acc_data_parsed,acc_context)) %>%
    setnames(c("id_acc","id","anno","Nmet","N","value")) %>%
    .[id%in%acc_feature] %>% .[N>=min_gpc]
  
  # Load RNA data
  sce <- load_SingleCellExperiment(io$rna.sce)[rna_gene,]
  # sce <- sce[rowData(sce)$symbol == rna_gene,]
  rna_dt <- logcounts(sce) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
    melt(id.vars = "id_rna", value.name = "value", variable.name = "id")
  rna_dt$id <- rna_gene
    # melt(id.vars = "id_rna", value.name = "value", variable.name = "ens_id") %>%
    # merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","id"))
  
  
  return(list("met"=met_dt, "acc"=acc_dt, "rna"=rna_dt))
}