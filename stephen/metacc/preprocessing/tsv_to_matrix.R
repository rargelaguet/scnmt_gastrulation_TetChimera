library(data.table)
library(purrr)
library(Matrix)



# script to process cpg/gpc level data into mean methylation per locus using a set of features


source(here::here("settings.R"))


opts$min_met_cov <- 1
opts$min_acc_cov <- 1

dir(io$features.dir)

opts$anno <- "H3K27ac_distal_E7.5_Mes_intersect12" #"multiome_peaks"


anno <- paste0(io$features.dir, "/", opts$anno, ".bed.gz") %>% 
  map(fread) %>%
  rbindlist() %>% 
  setnames(c("chr", "start", "end", "strand", "id", "anno"))

# convert chromosome column
if (!anno[1, chr %like% "chr"]) anno[, chr := paste0("chr", chr)]

# convert id to chr:start-end
anno[, id := paste0(chr, ":", start, "-", end)]

setkey(anno, chr, start, end)

metadata <- fread(io$metadata)

met_files <- metadata[pass_metQC == TRUE, paste0(io$met_data_raw, "/", cg_files)]
acc_files <- metadata[pass_accQC == TRUE, paste0(io$acc_data_raw, "/", gc_files)]


process_file <- function(file, anno){
  cell <- gsub(".tsv.gz", "", basename(file))
  
  dt <- fread(file)
  
  # convert rate to %
  if (dt[, max(rate)==1]) dt[, rate := rate * 100]
  
  # run overlap and compute mean methylation per position
  setkey(dt, chr, start, end) %>% 
    foverlaps(anno, nomatch = 0L) %>% 
    .[, .(cell = cell, rate = round(mean(rate)), .N), .(id, anno)]
}


met <- map(met_files, process_file, anno = anno) %>% rbindlist()
met

acc <- map(acc_files, process_file, anno = anno) %>% rbindlist()
acc

outfiles <- paste0(
  c(io$met_data_parsed, io$acc_data_parsed),
  "/",
  opts$anno,
  ".tsv.gz"
)
outfiles
walk2(list(met, acc), outfiles, fwrite_tsv)


# convert_to_matrix <- function(x, anno){
#   anno <- anno[, i := .I]
#   cells <- x[, .(cell = unique(cell))][, j := .I]
#   
#   positions <- x[N >= opts$min_met_cov, .(id, cell, rate)] %>% 
#     merge(anno[, .(id, i)], by = "id") %>% 
#     merge(cells, by = "cell")
#   
#   ids <- positions[, .(i, id)] %>% unique() %>% .[order(i)]
#   ids
#   
#   mat <- sparseMatrix(i = positions$i, j = positions$j, x = positions$rate)
#   
#   mat
#   
#   nrow(mat)
#   nrow(ids)
#   nrow(anno)
#   rownames(mat) <- ids$id
#   colnames(mat) <- cells$cell
# }
# 
# met_mat <- convert_to_matrix(met, anno)




