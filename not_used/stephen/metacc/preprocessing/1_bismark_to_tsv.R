library(data.table)
library(purrr)



# script to process bismark output files
# takes as input the bismark output folder(s)
# merges the read1 and read2 files, filters non-binary sites and outputs to tsv
# processes CpG methylation and GpC accessibility files separately 

# TODO - read sierra dirs from a file

source(here::here("settings.R"))

# in/out



# io$indir <- c("/bi/sequencing/Sample_5294_scNMT_plate1/Lane_7466_scNMT_plate1/Aligned/",
#               "/bi/sequencing/Sample_5295_scNMT_plate2/Lane_7467_scNMT_plate2/Aligned/",
#               "/bi/sequencing/Sample_5296_scNMT_plate3/Lane_7468_scNMT_plate3/Aligned/",
#               "/bi/sequencing/Sample_5297_scNMT_plate4/Lane_7469_scNMT_plate4/Aligned/",
#               "/bi/sequencing/Sample_5298_scNMT_plate5/Lane_7470_scNMT_plate5/Aligned/",
#               "/bi/sequencing/Sample_5299_scNMT_plate6/Lane_7471_scNMT_plate6/Aligned/")

io$indir <- c(
  "/bi/sequencing/Sample_5515_scNMT-1/Lane_7816_scNMT-1/Aligned/",
  "/bi/sequencing/Sample_5516_scNMT-2/Lane_7817_scNMT-2/Aligned/",
  "/bi/sequencing/Sample_5538_scNMT_CRI/Lane_7879_scNMT_CRI/Aligned/",
  "/bi/sequencing/Sample_5561_E8.5_oct20_plate4/Lane_7932_E8.5_oct20_plate4/Aligned/",
  "/bi/sequencing/Sample_5562_E8.5_oct20_plate5/Lane_7933_E8.5_oct20_plate5/Aligned/",
  "/bi/sequencing/Sample_5563_E8.5_oct20_plate6/Lane_7934_E8.5_oct20_plate6/Aligned/",
  "/bi/sequencing/Sample_5564_E8.5_oct20_plate7/Lane_7935_E8.5_oct20_plate7/Aligned/",
  "/bi/sequencing/Sample_5565_E8.5_oct20_plate8/Lane_7936_E8.5_oct20_plate8/Aligned/"
)




io$qc_out <- file.path(io$basedir, "metacc/qc/")

opts$filter_non_binary <- TRUE

# create folders
walk(list(io$met_data_raw, io$acc_data_raw, io$qc_out), dir.create, recursive = TRUE)


# find cpg and gpc files
cg <- dir(io$indir, pattern = "CpG.cov", full = TRUE)
gc <- dir(io$indir, pattern = "GpC.cov", full = TRUE)
in_files <- list(CpG = cg, GpC = gc)

cell_names <- map(in_files, ~gsub("_R[1-4].*", "", basename(.x)))


valid_chrs <- paste0("chr", c(1:19, "X", "Y"))




print("processing CpG methylation files....")

walk(unique(cell_names$CpG), ~{
  cell <- .x
  files <- in_files$CpG[cell_names$CpG == cell]
  
  print(paste("processing", .x))
  
  dt <- map(files, fread, select = c(1:2, 5:6)) %>% 
    rbindlist() %>% 
    .[, .(metreads = sum(V5), nonmetreads = sum(V6)), .(V1, V2)] %>% 
    .[, .(
      chr = V1,
      start = V2,
      end = V2,
      rate = metreads / (metreads + nonmetreads)
    )]
  
  # rename chromosomes
  if (dt[1, !chr %like% "chr"]) dt[, chr := paste0("chr", chr)]
  
  # filter non valid chrs
  dt <- dt[chr %in% valid_chrs]
  
  # filter non binary
  if (opts$filter_non_binary) dt <- dt[rate == 0 | rate == 1]
  
  # sort
  setkey(dt, chr, start)
  
  
  
  
  print("saving tsv...")
  
  # save as bigwig
  outfile <- paste0(io$met_data_raw, "/", .x, "_CpG.tsv.gz")
  fwrite_tsv(dt, outfile)
  
  print("saving QC...")
  
  # save QC info
  qc <- dt[, .(sample = cell, Ncpg = .N, mean_met = mean(rate), tsv = outfile, date_processed = Sys.Date())]
  qc_file <- paste0(io$qc_out, "/CpG_QC_", .x, ".tsv")
  fwrite(qc, qc_file)
})




print("processing GpC accessibility files....")

walk(unique(cell_names$GpC), ~{
  cell <- .x
  files <- in_files$GpC[cell_names$GpC == cell]
  
  print(paste("processing", .x))
  
  dt <- map(files, fread, select = c(1:2, 5:6)) %>% 
    rbindlist() %>% 
    .[, .(metreads = sum(V5), nonmetreads = sum(V6)), .(V1, V2)] %>% 
    .[, .(
      chr = V1,
      start = V2,
      end = V2,
      rate = metreads / (metreads + nonmetreads)
    )]
  
  # rename chromosomes
  if (dt[1, !chr %like% "chr"]) dt[, chr := paste0("chr", chr)]
  
  # filter non valid chrs
  dt <- dt[chr %in% valid_chrs]
  
  # filter non binary
  if (opts$filter_non_binary) dt <- dt[rate == 0 | rate == 1]
  
  # sort
  setkey(dt, chr, start)
  
  
  print("saving tsv...")
  
  # save as bigwig
  outfile <- paste0(io$acc_data_raw, "/", .x, "_GpC.tsv.gz")
  fwrite_tsv(dt, outfile)
  
  print("saving QC...")
  
  # save QC info
  qc <- dt[, .(sample = cell, Ngpc = .N, mean_acc = mean(rate), tsv = outfile, date_processed = Sys.Date())]
  qc_file <- paste0(io$qc_out, "/GpC_QC_", .x, ".tsv")
  fwrite(qc, qc_file)
})






