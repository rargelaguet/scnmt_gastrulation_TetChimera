library(data.table)
library(purrr)
library(Rsubread)


io <- list()
io$indirs <- c("/bi/sequencing/Sample_5223_scRNA-seq_1/Lane_7352_scRNA-seq_1/Aligned/",
               "/bi/sequencing/Sample_5224_scRNA-seq_2/Lane_7353_scRNA-seq_2/Aligned/",
               "/bi/sequencing/Sample_5225_scRNA-seq_3/Lane_7354_scRNA-seq_3/Aligned/",
               "/bi/sequencing/Sample_5226_scRNA-seq_4/Lane_7355_scRNA-seq_4/Aligned/",
               "/bi/sequencing/Sample_5227_scRNA-seq_5/Lane_7356_scRNA-seq_5/Aligned/",
               "/bi/sequencing/Sample_5228_scRNA-seq_6/Lane_7357_scRNA-seq_6/Aligned/")

#io$gtf <- "/bi/scratch/Stephen_Clark/gastrulation_data/features/genes/Mus_musculus.GRCm38.87.gtf"
io$gtf <- "/bi/scratch/Stephen_Clark/tet_chimera_nmtseq/features/gtf/Mus_musculus.GRCm38.100.gtf.gz"


io$outdir <- "/bi/scratch/Stephen_Clark/tet_chimera_nmtseq/rna/counts"

dir.create(io$outdir, recursive = TRUE)

infiles <- dir(io$indirs, full = TRUE, pattern = ".bam$") %>% unlist()
samples <- gsub(".bam", "", basename(infiles))

fc <- featureCounts(files = infiles,
                    annot.ext = io$gtf,
                    isGTFAnnotationFile = TRUE,
                    countMultiMappingReads = FALSE,
                    isPairedEnd = TRUE,
                    nthreads = 8)


counts <- as.data.frame(fc$counts) %>%
  setDT(keep.rownames = "ens_id") %>%
  setnames(colnames(.), gsub("\\.", "_", colnames(.)))

meta <- data.table(id_rna = colnames(counts)[2:ncol(counts)]) %>%
  .[, plate := strsplit(id_rna, "_") %>% map_chr(1)] %>%
  .[, well := strsplit(id_rna, "_") %>% map_chr(2)]
  

fwrite(counts, paste0(io$outdir, "/counts_GRCm38_100.tsv.gz"), sep = "\t", na = "NA")
fwrite(meta, paste0(io$outdir, "/cellnames.tsv.gz"), sep = "\t", na = "NA")
