library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)

source(here::here("settings.R"))

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






