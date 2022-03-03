library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)

source(here::here("settings.R"))

io$features <- "/bi/scratch/Stephen_Clark/multiome/raw/processed/atac/archR/peakSet/peakSet.rds"

infiles <- dir(
  "/bi/scratch/Stephen_Clark/tet_chimera_nmtseq/met/raw/bigwigs/",
  pattern = ".bw$",
  full = TRUE
)

features <- readRDS(io$features)
features$anno <- names(features)
selection <- BigWigSelection(ranges = features) # should be able to use this to import only the ranges we want. but its slower than importing and runing overlap

.x=files[1]
gr <- import.bw(.x)

ol <- findOverlaps(gr, features)

