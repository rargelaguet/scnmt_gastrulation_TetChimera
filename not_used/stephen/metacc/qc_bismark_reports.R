library(data.table)
library(purrr)

source(here::here("settings.R"))


# TODO: create plate metadata file with sierra locations


io$bismark_summary <- c(
  "/bi/sequencing/Sample_5515_scNMT-1/Lane_7816_scNMT-1/Aligned/plate1_L001_bismark_summary_report.txt",
  "/bi/sequencing/Sample_5516_scNMT-2/Lane_7817_scNMT-2/Aligned/plate2_L002_bismark_summary_report.txt"
)

io$oufile <- paste0(io$basedir, "metacc/bismark_qc.tsv.gz")

summary <- map(io$bismark_summary, fread) %>% 
  map2(basename(dirname(dirname(io$bismark_summary))), ~.x[, folder := .y]) %>% 
  rbindlist() %>% 
  setnames(colnames(.), make.names(colnames(.)))

summary
