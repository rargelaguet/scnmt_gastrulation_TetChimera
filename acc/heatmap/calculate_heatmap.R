# Define settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera/acc/heatmap/load_settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation_TetChimera/acc/heatmap/load_settings.R")
} else {
  stop("Computer not recognised")
}

# Load data and summarise
dt <- lapply(opts$annos, function(n) {
  print(n)
  file = sprintf("%s/%s.tsv.gz",io$acc_data_parsed,n)
  if (file.exists(file)) {
    fread(file, showProgress = F, header = F,
      select = c("V1"="character","V2"="character","V3"="factor","V4"="integer","V5"="integer","V6"="integer")
      ) %>% setnames(c("id_acc","id","anno","Nmet","Ntotal","rate")) %>%
      .[,.(rate=100*(sum(Nmet)/sum(Ntotal))),by = c("id_acc","anno")]
  }
}) %>% rbindlist

# Save
fwrite(dt, paste0(io$outdir,"/heatmap_acc.tsv.gz"), sep="\t")
