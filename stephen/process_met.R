library(data.table)
library(purrr)

source(here::here("functions.R"))

indir <- c("/bi/sequencing/Sample_5104_scNMT_Tet_chimera_1/Lane_7147_scNMT_Tet_chimera_1/Aligned/",
           "/bi/sequencing/Sample_5105_scNMT_Tet_chimera_2/Lane_7148_scNMT_Tet_chimera_2/Aligned/",
           "/bi/sequencing/Sample_5106_scNMT_Tet_chimera_3/Lane_7149_scNMT_Tet_chimera_3/Aligned/",
           "/bi/sequencing/Sample_5107_scNMT_Tet_chimera_4/Lane_7150_scNMT_Tet_chimera_4/Aligned/")

cg_files <- dir(indir, pattern = "bismark.cov.gz", f = T) %>% unlist()

names <- basename(cg_files) %>% gsub("_L00.*", "", .)

files <- data.table(file = cg_files, name = names, read = rep(c("file1", "file2"), length(cg_files)/2)) %>%
  dcast(name ~ read, value.var = "file" )

files=files[sample(1:.N, 20)]

met <- pmap(files, read_covs)
