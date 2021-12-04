library(data.table)
library(purrr)
library(furrr)

# TODO - read sierra dirs from a file

source(here::here("metacc/preprocessing/BS_preprocessing_functions.R"))

# in/out

io <- list()
io$data_dir <- "/bi/scratch/Stephen_Clark/tet_chimera_nmtseq"

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



io$cg_outdir <- file.path(io$data_dir, "met/raw/oct20")
io$gc_outdir <- file.path(io$data_dir, "acc/raw/oct20")

io$qc_out <- file.path(io$data_dir, "metacc/qc/")
    
# create folders
create_dir(io$cg_outdir)
create_dir(io$gc_outdir)
create_dir(io$qc_out)

cg <- dir(io$indir, pattern = "CpG.cov", full = TRUE)
gc <- dir(io$indir, pattern = "GpC.cov", full = TRUE)




in_files <- list(CpG = cg, GpC = gc)

out_dt <- data.table(
  in_files = in_files, 
  outdir = c(io$cg_outdir, io$gc_outdir),
  suffix = c("_CpG.tsv.gz", "_GpC.tsv.gz")
  )

out_files <- pmap(out_dt, function(in_files, outdir, suffix){
  paste0(outdir, "/", gsub("_R[1-4].*", "", basename(in_files)), suffix)
}) %>%
  set_names(names(in_files))
  





qc <- process_bismark_files(
  in_files = in_files,
  out_files = out_files,
  columns = c(1:2, 5:6), 
  binarise = TRUE,
  parallel = FALSE, 
  sort = TRUE,
  gzip = FALSE
)



  # save QC file

qc_file <- paste0(io$qc_out, "/", Sys.Date(), ".tsv")
i <- 1
while (file.exists(qc_file)) {
    qc_file <- gsub("[0-9].tsv$|.tsv$", paste0(i, ".tsv"), qc_file)
    i <- i + 1
}

fwrite_tsv(qc, qc_file)
