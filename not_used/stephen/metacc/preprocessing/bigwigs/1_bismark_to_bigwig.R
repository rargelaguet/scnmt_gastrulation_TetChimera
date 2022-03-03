library(data.table)
library(purrr)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)


# script to process bismark output files
# takes as input the bismark output folder(s)
# merges the read1 and read2 files, filters non-binary sites and outputs to bigwig
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
  "/bi/sequencing/Sample_5516_scNMT-2/Lane_7817_scNMT-2/Aligned/"
)



io$cg_outdir <- file.path(io$basedir, "met/raw/bigwigs")
io$gc_outdir <- file.path(io$basedir, "acc/raw/bigwigs")

io$qc_out <- file.path(io$basedir, "metacc/qc/")

opts$filter_non_binary <- TRUE

# create folders
walk(list(io$cg_outdir, io$gc_outdir, io$qc_out), dir.create, recursive = TRUE)


# find cpg and gpc files
cg <- dir(io$indir, pattern = "CpG.cov", full = TRUE)
gc <- dir(io$indir, pattern = "GpC.cov", full = TRUE)
in_files <- list(CpG = cg, GpC = gc)

cell_names <- map(in_files, ~gsub("_R[1-4].*", "", basename(.x)))

out_files <- pmap(out_dt, function(in_files, outdir, suffix){
  paste0(outdir, "/", gsub("_R[1-4].*", "", basename(in_files)), suffix)
}) %>%
  set_names(names(in_files))

valid_chrs <- paste0("chr", c(1:19, "X", "Y"))
seqs <- seqinfo(BSgenome.Mmusculus.UCSC.mm10) %>% .[valid_chrs]



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
      score = metreads / (metreads + nonmetreads)
    )]
  
  # rename chromosomes
  if (dt[1, !chr %like% "chr"]) dt[, chr := paste0("chr", chr)]
  
  # filter non valid chrs
  dt <- dt[chr %in% valid_chrs]
  
  # filter non binary
  if (opts$filter_non_binary) dt <- dt[score == 0 | score == 1]
  
  # sort
  setkey(dt, chr, start)
  
  # convert to Granges
  gr <- makeGRangesFromDataFrame(
    dt, 
    keep.extra.columns = TRUE,
    seqinfo = seqs
  )
  
  print("saving bigwig...")
  
  # save as bigwig
  outfile <- paste0(io$cg_outdir, "/", .x, "_CpG.bw")
  export.bw(gr, outfile)
  
  print("saving QC...")
  
  # save QC info
  qc <- dt[, .(sample = cell, Ncpg = .N, mean_met = mean(score), bigwig = outfile, date_processed = Sys.Date())]
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
      score = metreads / (metreads + nonmetreads)
    )]
  
  # rename chromosomes
  if (dt[1, !chr %like% "chr"]) dt[, chr := paste0("chr", chr)]
  
  # filter non valid chrs
  dt <- dt[chr %in% valid_chrs]
  
  # filter non binary
  if (opts$filter_non_binary) dt <- dt[score == 0 | score == 1]
  
  # sort
  setkey(dt, chr, start)
  
  # convert to Granges
  gr <- makeGRangesFromDataFrame(
    dt, 
    keep.extra.columns = TRUE,
    seqinfo = seqs
  )
  
  print("saving bigwig...")
  
  # save as bigwig
  outfile <- paste0(io$gc_outdir, "/", .x, "_GpC.bw")
  export.bw(gr, outfile)
  
  print("saving QC...")
  
  # save QC info
  qc <- dt[, .(sample = cell, Ngpc = .N, mean_acc = mean(score), bigwig = outfile, date_processed = Sys.Date())]
  qc_file <- paste0(io$qc_out, "/GpC_QC_", .x, ".tsv")
  fwrite(qc, qc_file)
})






