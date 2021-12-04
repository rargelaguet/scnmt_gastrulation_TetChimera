library(data.table)
library(purrr)
library(Rsubread)


io <- list()
io$indirs <- c(
  "/bi/sequencing/Sample_5381_E7.5_tet_crispr_RNA/Lane_7633_E7.5_tet_crispr_RNA/Aligned/",
  "/bi/sequencing/Sample_5448_tet_chimera_rna_pool1/Lane_7710_tet_chimera_rna_pool1/Aligned/",
  "/bi/sequencing/Sample_5453_tet_chimera_rna_pool2/Lane_7726_tet_chimera_rna_pool2/Aligned/"
  )

#io$gtf <- "/bi/scratch/Stephen_Clark/gastrulation_data/features/genes/Mus_musculus.GRCm38.87.gtf"
io$gtf <- "/bi/scratch/Stephen_Clark/tet_chimera_nmtseq/features/gtf/Mus_musculus.GRCm38.100.gtf.gz"
io$outdir <- "/bi/scratch/Stephen_Clark/tet_chimera_nmtseq/rna/counts/new2021"


opts <- list()
opts$cores <- 32
opts$paired <- FALSE


dir.create(io$outdir, recursive = TRUE)

infiles <- dir(io$indirs, full = TRUE, pattern = ".bam$") %>% unlist()


samples <- gsub(".bam", "", basename(infiles))


# select <- samples %like% "crispr"
# table(select)
# 
# infiles <- infiles[select]
# samples <- samples[select]

fc <- featureCounts(
  files                  = infiles,
  annot.ext              = io$gtf,
  isGTFAnnotationFile    = TRUE,
  countMultiMappingReads = FALSE,
  isPairedEnd            = opts$paired,
  nthreads               = opts$cores
  )


counts <- as.data.frame(fc$counts) %>%
  setDT(keep.rownames = "ens_id") %>%
  setnames(c("ens_id", samples))

# meta <- data.table(id_rna = colnames(counts)[2:ncol(counts)]) %>%
#   .[, plate := strsplit(id_rna, "_") %>% map_chr(1)] %>%
#   .[, well := strsplit(id_rna, "_") %>% map_chr(2)]
  

fwrite(counts, paste0(io$outdir, "/counts.tsv.gz"), sep = "\t", na = "NA", quote = FALSE)
# fwrite(meta, paste0(io$outdir, "/meta.tsv.gz"), sep = "\t", na = "NA", quote = FALSE)
