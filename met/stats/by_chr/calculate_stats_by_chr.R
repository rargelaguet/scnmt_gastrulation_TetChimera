#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera//settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation_TetChimera//settings.R")
} else {
  stop("Computer not recognised")
}

# Define I/O
io$outdir <- paste0(io$basedir,"/met/results/stats")

# Define cells
opts$cells <- list.files(io$met_data_raw, pattern = "*.tsv.gz") %>% 
  stringr::str_replace_all(".tsv.gz","")

#####################
## Calculate stats ##
#####################

stats <- data.table(expand.grid(opts$cells,opts$chr)) %>% 
  setnames(c("id_met","chr")) %>%
  .[,c("nCG","mean"):=as.numeric(NA)]

for (i in opts$cells) {
  if (file.exists(sprintf("%s/%s.tsv.gz",io$met_data_raw,i))) {
    print(i)

    # Load sample methylation data
    data <- fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,i), sep=" ", verbose=F, showProgress=F) %>%
      setnames(c("chr","pos","rate")) %>%
      .[,chr:=paste0("chr",chr)]

    # Compute methylation statistics per chromosome
    for (j in opts$chr) {
      data_j <- data[chr==j]
      stats[id_met==i & chr==j, c("nCG","mean"):=list(nrow(data_j),round(mean(data_j$rate)*100,2))]  
    }

  } else {
    print(sprintf("Sample %s not found",i))
  }
}

# Save
fwrite(stats, paste0(io$outdir,"/stats_per_chromosome.txt.gz"), sep="\t", na="NA")
