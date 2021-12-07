suppressMessages(library(argparse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--indir',  type="character",              help='Input directory')
p$add_argument('--outdir',  type="character",              help='Output directory')
p$add_argument('--context',  type="character",              help='cg/CG or gc/GC')
p$add_argument('--test', action="store_true",             help='Test mode? subset number of cells')


# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args <- list()
# args$indir <- file.path(io$basedir,"processed/met/cpg_level")
# args$outdir <- file.path(io$basedir,"processed/met/feature_level")
# args$featuresdir  <- file.path(io$basedir,"features/genomic_contexts")
# args$metadata <- file.path(io$basedir,"results/met/qc/sample_metadata_after_met_qc.txt.gz")
# args$context <- "CG"
# args$annos <- c("prom_2000_2000")
## END TEST ##

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))




##############
## ORIGINAL ##
##############

io$outdir <- "/bi/scratch/Stephen_Clark/tet_chimera_nmtseq/met/pseudobulk"
dir.create(io$outdir, recursive = TRUE)

#metadata <- fread(io$metadata)

infiles <- dir(
  "/bi/scratch/Stephen_Clark/tet_chimera_nmtseq/met/raw/bigwigs/",
  pattern = ".bw$",
  full = TRUE
) %>% 
  data.table(file = ., plate = strsplit(., "_") %>% map_chr(6)) %>% 
  split(by = "plate")

valid_chrs <- paste0("chr", c(1:19, "X", "Y"))
seqs <- seqinfo(BSgenome.Mmusculus.UCSC.mm10) %>% .[valid_chrs]


walk(names(infiles), ~{
  files <- infiles[[.x]][, file]
  
  gr <- map(files, import.bw) %>% 
    map(as.data.table) %>% 
    map(~.x[, .(
      id = paste(seqnames, start, end, sep = "_"),
      score
    )]) %>% 
    map2(basename(files), ~setnames(.x, "score", .y)) %>% 
    purrr::reduce(merge, by = c("id"), all = TRUE) %>% 
    .[, .(id = id, score = apply(.SD, 1, mean, na.rm = TRUE)), .SDcol = colnames(.)[2:ncol(.)]] %>% 
    tidyr::separate("id", c("chr", "start", "end")) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = seqs)
  
  outfile <- paste0(io$outdir, "/", .x, ".bw")
  export.bw(gr, outfile)
})






