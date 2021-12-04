library(data.table)
library(purrr)




# script to process cpg/gpc level data into mean methylation per locus using a set of features


source(here::here("settings.R"))

# dir(io$features.dir)


opts$min_met_cov <- 1
opts$min_acc_cov <- 1


opts$anno <- c(
  "100000bp_bins_50000_step",
  "10000bp_bins_5000_step"
)

print(paste("features selected:", opts$anno))


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

#met_files <- metadata[pass_metQC == TRUE, paste0(io$met_data_raw, "/", cg_files)]
#acc_files <- metadata[pass_accQC == TRUE, paste0(io$acc_data_raw, "/", gc_files)]
met_files <- metadata[pass_metQC == TRUE, cg_files]
acc_files <- metadata[pass_accQC == TRUE, gc_files]



process_file <- function(file, anno){
  
  cell <- gsub(".tsv.gz", "", basename(file))
  
  print(paste("processing", cell))
  
  dt <- fread(file)
  
  # convert rate to %
  if (dt[, max(rate)==1]) dt[, rate := rate * 100]
  
  # run overlap and compute mean methylation per position
  setkey(dt, chr, start, end) %>% 
    foverlaps(anno, nomatch = 0L) %>% 
    .[, .(cell = cell, rate = round(mean(rate)), .N), .(id, anno)]
}


met <- map(met_files, process_file, anno = anno) %>% 
  rbindlist() %>% 
  split(by = "anno", keep.by = TRUE)
head(met)


# save
print("saving methylation file...")

metfiles <- paste0(
  io$met_data_parsed,
  "/",
  names(met),
  ".tsv.gz"
)

print(metfiles)
walk2(met, metfiles, fwrite_tsv)


acc <- map(acc_files, process_file, anno = anno) %>% 
  rbindlist()%>% 
  split(by = "anno", keep.by = TRUE)
head(acc)


print("saving accessibility file...")


accfiles <- paste0(
  io$acc_data_parsed,
  "/",
  names(acc),
  ".tsv.gz"
)



print(accfiles)
walk2(acc, accfiles, fwrite_tsv)






