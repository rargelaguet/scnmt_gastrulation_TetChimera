source(here::here("settings.R"))


io$de_genes <- c(
  paste0(io$basedir, "/other/10x_DE/E7.5/"),
  paste0(io$basedir, "/other/10x_DE/E8.5/")
  
)

io$outfile <- gsub("\\..*", "_DEgenes.tsv.gz", io$features.tsv)

# how many bp up/downstream of gene
opts$extend_by <- list(
  upstream   = 5e4,
  downstream = 5e4
)



# load gene metadata and extend start/end positions
gene_metadata <- fread(io$gene_metadata) %>% 
  .[strand == "+", start := start - opts$extend_by$upstream] %>% 
  .[strand == "+", end   := end   + opts$extend_by$downstream] %>% 
  .[strand == "-", start := start - opts$extend_by$downstream] %>% 
  .[strand == "-", end   := end   + opts$extend_by$upstream] %>% 
  .[start < 0, start := 0] # if extension goes beyond start of chromosome

# load DE genes from 10x tet chimera data
de_genes <- dir(io$de_genes, full = TRUE, pattern = ".txt.gz$") %>% 
  set_names(gsub(".txt.gz$", "", basename(.))) %>% 
  map(fread) %>% 
  map2(names(.), ~.x[, celltype := .y]) %>% 
  rbindlist(use.names = TRUE, fill = TRUE) %>% 
  .[sig == TRUE]


# load features
features <- fread(io$features.tsv) %>% 
  setnames(c("chr", "start", "end", "strand", "id", "anno"))

if (!features[1, chr %like% "chr"]) features[, chr := paste0("chr", chr)]

head(features)


setkey(features, chr, start, end)


# merge DE genes with gene metadata then run overlap with features

overlap <- de_genes[, .(ens_id, updown = sign(logFC), celltype)] %>% 
  merge(gene_metadata, by = "ens_id") %>% 
  setkey(chr, start, end) %>% 
  foverlaps(features, nomatch = 0L)


fwrite_tsv(overlap, io$outfile)



