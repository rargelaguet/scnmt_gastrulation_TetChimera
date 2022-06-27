library(data.table)
library(purrr)

# this script pulls the bismark cytosine context summary files from the bi sqeuencing folders
# then filters to include only WCH trinucleotides
# then computes a conversion rate estimate for each cell

outfile <- "/bi/scratch/Stephen_Clark/tet_chimera_nmtseq/metacc/qc/conversion/conversion.tsv.gz"
dir.create(dirname(outfile))

seqdirs <- c(
  "/bi/sequencing/Sample_5515_scNMT-1/Lane_7816_scNMT-1/Aligned/",
  "/bi/sequencing/Sample_5516_scNMT-2/Lane_7817_scNMT-2/Aligned/",
  "/bi/sequencing/Sample_5538_scNMT_CRI/Lane_7879_scNMT_CRI/Aligned/",
  "/bi/sequencing/Sample_5561_E8.5_oct20_plate4/Lane_7932_E8.5_oct20_plate4/Aligned/",
  "/bi/sequencing/Sample_5562_E8.5_oct20_plate5/Lane_7933_E8.5_oct20_plate5/Aligned/",
  "/bi/sequencing/Sample_5563_E8.5_oct20_plate6/Lane_7934_E8.5_oct20_plate6/Aligned/",
  "/bi/sequencing/Sample_5564_E8.5_oct20_plate7/Lane_7935_E8.5_oct20_plate7/Aligned/",
  "/bi/sequencing/Sample_5565_E8.5_oct20_plate8/Lane_7936_E8.5_oct20_plate8/Aligned/"
)


files <- dir(seqdirs, pattern = "cytosine_context_summary.txt", full = TRUE)
cells <- strsplit(basename(files), "_") %>% map_chr(~paste(.x[1], .x[2], .x[3], .x[4], .x[5], sep = "_"))

stats <- map2(files, cells, ~fread(.x)[, cell := .y]) %>% 
  rbindlist()
stats

stats[, unique(`full context`)]

stats[, first := substr(`full context`, 1, 1)]
stats[, third := substr(`full context`, 3, 3)]

filtered <- stats[first %in% c("A", "T") & !third %in% "G"]

conversion <- filtered[, .(met = sum(`count methylated`), unmet = sum(`count unmethylated`)), cell] %>% 
  .[, conversion := unmet / (met + unmet)]

fwrite(conversion, outfile, sep = "\t", quote = FALSE, na = "NA")
