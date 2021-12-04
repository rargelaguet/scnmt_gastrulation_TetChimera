source(here::here("settings.R"))


opts$min_cg_cov  <- 1e5
opts$min_gc_rate <- 0.1

io$outfile <- paste0(io$basedir, "/metacc/qc/merged/metacc_QC.tsv.gz")


extract_name <- function(x){
  sp <- strsplit(x, "_") %>% unlist()
  stage <- sp[sp %like% "E[1-9]\\.5"]
  batch <- sp[sp %like% "tet|chimera|oct20"] %>% paste(collapse = "_")
  plate <- sp[sp %like% "plate|Plate"]
  well <- sp[sp %like% "[A-H][01-12]"]
  paste(stage, batch, plate, well, sep = "_")
}

# extract_name <- function(x){
#   gsub("L00[1-9]", "", x) %>% 
#     gsub("[A-H][01-12].*", "", .)
# }

##### not working yet

cpg_qc <- dir(
  paste0(io$basedir, "/metacc/qc"), 
  pattern = "CpG", 
  full = TRUE
) %>% 
  map(fread) %>% 
  rbindlist() %>% 
  .[, cell := map_chr(sample, extract_name)] %>% 
  .[, pass_metQC := ifelse(Ncpg >= opts$min_cg_cov, TRUE, FALSE)] %>% 
  setnames("tsv", "cg_files")

gpc_qc <- dir(
  paste0(io$basedir, "/metacc/qc"), 
  pattern = "GpC", 
  full = TRUE
) %>% 
  map(fread) %>% 
  rbindlist()%>% 
  .[, cell := map_chr(sample, extract_name)] %>% 
  .[, pass_accQC := ifelse(mean_acc >= opts$min_gc_rate, TRUE, FALSE)] %>% 
  setnames("tsv", "gc_files")

combined_qc <- merge(cpg_qc, gpc_qc, by = "cell") %>% 
  .[, .(
    cell,
    Ncpg,
    mean_met,
    pass_metQC,
    Ngpc,
    mean_acc,
    pass_accQC = as.logical(pass_metQC * pass_accQC), # only pass accQC if also meets min cpgs
    cg_files,
    gc_files
  )]

fwrite_tsv(combined_qc, io$outfile)
#combined_qc=fread(io$outfile)

met_stats <- cpg_qc[, .(
  id_met = basename(cg_files) %>% gsub(".tsv.gz", "", .),
  Ncpg,
  mean = mean_met
)]

acc_stats <- gpc_qc[, .(
  id_acc = basename(gc_files) %>% gsub(".tsv.gz", "", .),
  Ngpc,
  mean = mean_acc
)]

fwrite_tsv(met_stats, io$met.stats)
fwrite_tsv(acc_stats, io$acc.stats)

sample_metadata <- fread(io$metadata)

cols <- colnames(sample_metadata) %>% 
  .[. %in% c("pass_accQC", "pass_metQC", "cg_files", "gc_files")]

sample_metadata[, c(cols) := NULL]


sample_metadata <- merge(
  sample_metadata, 
  combined_qc[, .(cell, pass_metQC, pass_accQC, cg_files, gc_files)], 
  by = "cell", 
  all.x = TRUE
)

fwrite_tsv(sample_metadata, io$metadata)


  # library(ggplot2)
  # library(cowplot)
  # ggplot(combined_qc, aes(Ncpg)) + 
  #   geom_histogram() +
  #   theme_cowplot()
