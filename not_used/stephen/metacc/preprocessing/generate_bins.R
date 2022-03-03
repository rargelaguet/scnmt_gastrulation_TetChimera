library(data.table)
library(purrr)




# script to process cpg/gpc level data into mean methylation per locus using a set of features



source(here::here("settings.R"))

io$chr_sizes <- "https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"

io$outdir <- io$features.dir


opts$binsize <- 1e5
opts$step    <- 5e4


feature_name <- paste0(
  format(opts$binsize, scientific = FALSE), 
  "bp_bins_",
  format(opts$step, scientific = FALSE),
  "_step"
)
print(feature_name)

sizes <- fread(io$chr_sizes) %>% 
  setnames(c("chr", "length")) %>% 
  .[!chr %like% "chrUn|random"]


bins <- sizes[, .(start = seq(0, length, by = opts$step)), chr] %>% 
  .[, end := start + opts$binsize] %>% 
  .[, c("strand", "id", "anno") := .("*", paste0(feature_name, "_", .I), feature_name)]

head(bins)

outfile <- paste0(io$outdir, "/", feature_name,".bed.gz")
print(paste("saving", outfile))

fwrite_tsv(bins, outfile, col.names = FALSE)

