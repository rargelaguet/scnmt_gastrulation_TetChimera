library(data.table)
library(purrr)
library(Rsubread)


io <- list()
io$indirs <- c("/bi/sequencing/Sample_4165_NMT_embryo_Transcriptome_Pool_1/Lane_5753_NMT_embryo_Transcriptome_Pool_1/Aligned/",
               "/bi/sequencing/Sample_4166_NMT_embryo_Transcriptome_Pool_2/Lane_5754_NMT_embryo_Transcriptome_Pool_2/Aligned/",
               "/bi/sequencing/Sample_4167_NMT_embryo_Transcriptome_Pool_3/Lane_5763_NMT_embryo_Transcriptome_Pool_3/Aligned/",
               "/bi/sequencing/Sample_4363_E4.5_E5.5_RNA/Lane_5998_E4.5_E5.5_RNA/Aligned/",                  
               "/bi/sequencing/Sample_4648_E7.5_endoderm_enriched_RNA/Lane_6442_E7.5_endoderm_enriched_RNA/Aligned/",     
               "/bi/sequencing/Sample_4682_E7.5_endoderm_enriched_RNA_p1p2/Lane_6504_E7.5_endoderm_enriched_RNA_p1p2/Aligned/",
               "/bi/sequencing/Sample_4738_embryos_dec18_PS_VE_RNA_pool1/Lane_6565_embryos_dec18_PS_VE_RNA_pool1/Aligned/",  
               "/bi/sequencing/Sample_4741_embryos_dec18_PS_VE_RNA_pool2/Lane_6585_embryos_dec18_PS_VE_RNA_pool2/Aligned/",  
               "/bi/sequencing/Sample_4753_embryos_dec18_PS_VE_RNA_pool2/Lane_6604_embryos_dec18_PS_VE_RNA_pool2/Aligned/")

#io$gtf <- "/bi/scratch/Stephen_Clark/gastrulation_data/features/genes/Mus_musculus.GRCm38.87.gtf" # Ricard's protein coding only annotation
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
  

fwrite(counts, paste0(io$outdir, "/gastrulation_counts_GRCm38_100.tsv.gz"), sep = "\t", na = "NA")
fwrite(meta, paste0(io$outdir, "/gastrulation_cellnames.tsv.gz"), sep = "\t", na = "NA")
