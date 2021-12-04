################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera//settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/scnmt_gastrulation_TetChimera//settings.R")
} else {
  stop("Computer not recognised")
}

# Output directory
io$outdir <- paste0(io$basedir,"/met/results/heatmap")

####################
## Define options ##
####################

# Define lineages
opts$lineages <- c(
  "hESC",
  "Zygote", 
  "2cell", 
  "4cell", 
  "8cell",
  "Morula", 
  "ICM", 
  "TE"
)

# Define genomic contexts for methylation
# opts$met.annos <- NULL
# opts$annos <- c("CGI","prom_2000_2000_cgi","genebody","LINE","L1","Alu","CR1")
if (is.null(opts$annos)) {
  opts$annos <- list.files(io$met_data_parsed, pattern=".tsv.gz") %>% gsub(".tsv.gz","",.)# %>% head(n=2)
}

#####################
## Update metadata ##
#####################

sample_metadata <- sample_metadata %>% 
  .[lineage%in%opts$lineage] %>%
  # .[sex%in%c("Female","Male")] %>%
  # .[pass_metQC==TRUE]
  .[!is.na(id_met)]
