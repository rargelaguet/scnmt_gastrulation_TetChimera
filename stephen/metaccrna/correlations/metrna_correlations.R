## ----echo=TRUE, include=FALSE--------------------------------------------------------------------------------------------
source(here::here("settings.R"))

library(SingleCellExperiment)

## ------------------------------------------------------------------------------------------------------------------------
io$outdir <- paste0(io$basedir,"/metaccrna/exploration/correlations")


## ------------------------------------------------------------------------------------------------------------------------
dir(io$features.dir)

# Define genomic contexts
opts$annos <- c(
  "multiome_peaks"
)

# how much upstream/downstream to include?
opts$extend_by <- list(
  upstream   = 5e4,
  downstream = 5e4
)

opts$split_by <- "class" # perform correlations separately split by this metadata variable

# Filtering parameters
opts$min.CpGs <- 1      # Minimum number of CpGs per feature
opts$min.cells <- 20    # Minimum number of observed cells per feature with at least opts$min.CpGs
opts$variance_filter <- 0.5  # keep the top x fraction of genes and features by variance

# Multiple testing correction options
opts$threshold_fdr  <- 0.1

# Correlation type options
opts$method <- "pearson"      # correlation type (see ?cor)
opts$weight <- FALSE          # weighted correlation (see ?wtd.cor) 

# Permutation test options
opts$permutation <- TRUE   # do permutation test?
opts$n_perms <- 10          # Number of random permutations

## ------------------------------------------------------------------------------------------------------------------------

# filter by QC
sample_metadata <- sample_metadata[pass_metQC == TRUE & pass_rnaQC == TRUE]


celltypes <- sample_metadata[, unique(celltype.mapped)]


## ------------------------------------------------------------------------------------------------------------------------

met.stats <- fread(io$met.stats) %>% .[,c("id_met","mean")] %>% .[,context:="CG"]



## ------------------------------------------------------------------------------------------------------------------------


# load gene metadata and extend start/end positions
gene_metadata <- fread(io$gene_metadata) %>% 
  .[strand == "+", start := start - opts$extend_by$upstream] %>% 
  .[strand == "+", end   := end   + opts$extend_by$downstream] %>% 
  .[strand == "-", start := start - opts$extend_by$downstream] %>% 
  .[strand == "-", end   := end   + opts$extend_by$upstream] %>% 
  .[start < 0, start := 0] # if extension goes beyond start of chromosome

# correct chromosome naming
if (!gene_metadata[1, chr %like% "chr"]) gene_metadata[, chr := paste0("chr", chr)]
setkey(gene_metadata, chr, start, end)

# load feature metadata 
feature_metadata <- paste0(io$features.dir, "/", opts$annos, ".bed.gz") %>% 
  map(fread) %>% 
  rbindlist() %>% 
  setnames(c("chr", "start", "end", "strand", "id", "anno"))

# correct chromosome naming
if (!feature_metadata[1, chr %like% "chr"]) feature_metadata[, chr := paste0("chr", chr)]

feature_metadata[,idx:=sprintf("%s:%s-%s",chr,start,end)]

# overlap with genes
setkey(feature_metadata, chr, start, end)

features_overlap <- foverlaps(gene_metadata, feature_metadata, nomatch = 0L) %>% 
  .[, .(id, anno, idx, ens_id, symbol)]

head(features_overlap)

## ----load_metdata-----------------------------------------------------------------------------

# read in methylation data, filter for QC, add cell info and gene info
met_dt <- paste0(io$met_data_parsed, "/", opts$annos, ".tsv.gz") %>% 
  map(fread) %>% 
  rbindlist() %>% 
  setnames("cell", "id_met") %>% 
  merge(sample_metadata[, .(id_met, cell, class, stage_lineage, markers)], by = "id_met") %>% 
  merge(features_overlap[, .(id = idx, anno, ens_id, gene = symbol)], by = c("id", "anno"), allow.cartesian = TRUE)




head(met_dt)


## ----load_rnadata-----------------------------------------------------------------------------

# load sce object
sce <- readRDS(io$rna.sce) 

# extract normalised expression counts, format for merging and filter QC
rna_dt <- sce@assays@data$logcounts %>% 
  as.data.table(keep.rownames = "ens_id") %>% 
  melt(id.vars = "ens_id", variable.name = "id_rna", value.name = "exp") %>% 
  merge(sample_metadata[, .(id_rna, cell)], by = "id_rna")

# Remove genes with little variability
rna_dt <- rna_dt[, var := var(expr), ens_id] %>% .[var>0.5] %>% .[,var:=NULL]

## ----join RNA to met-----------------------------------------------------------------------------


metrna <- merge(met_dt, rna_dt, by = c("ens_id", "cell"))
head(metrna)


## ----split data by phenotype--------------------------------------------------


metrna <- split(metrna, by = opts$split_by)

## ----filter by coverage then variance-----------------------------------------

# coverage
metrna <- map(metrna, ~{
  .x[N > opts$min.CpGs] %>% 
    .[, ncells := .N, .(id, anno, gene)] %>% 
    .[ncells >= opts$min.cells]
})

# filter by rna variance
metrna <- map(metrna, ~{
  genes <- .x[, .(cell, ens_id, exp)] %>% 
    unique() %>% 
    .[, .(var = var(exp)), .(ens_id)] %>% 
    setorder(-var) %>% 
    .[1: (.N * opts$variance_filter), ens_id]
  
  .x[ens_id %in% genes]
})

# filter by met variance
metrna <- map(metrna, ~{
  features <- .x[, .(cell, id, rate)] %>% 
    unique() %>% 
    .[, .(var = var(rate)), .(id)] %>% 
    setorder(-var) %>% 
    .[1: (.N * opts$variance_filter), id]
  
  .x[id %in% features]
})


## ----compute correlations-----------------------------------------------------

if (opts$weight){
  corfun <- function(exp, rate, N){
    x <- weights::wtd.cor(exp, rate, N)
    list(r = x[1], p = x[4])
  }
} else {
  corfun <- function(exp, rate, N){
    x <- cor.test(exp, rate)
    list(r = x$estimate, p = x$p.value)
  }
}

cors <- map(metrna, ~{
  .x[, corfun(exp, rate, N), .(id, ens_id, gene)] %>% 
    .[, c("padj_fdr", "padj_bonf") := list(p.adjust(p, method="fdr"), p.adjust(p, method="bonferroni"))] %>%
    .[, c("log_padj_fdr","log_padj_bonf") := list(-log10(padj_fdr), -log10(padj_bonf))] %>%
    .[, sig := padj_fdr <= opts$threshold_fdr] %>%  setorder(padj_fdr)
}) 

map(cors, ~.x[, sum(sig, na.rm = TRUE)])
  

cors <- map2(cors, names(cors), ~.x[, c(opts$split_by) := .y]) %>% 
  rbindlist()

outfile <- paste0(io$outdir, "/correlations.tsv.gz")
fwrite_tsv(cors, outfile)
