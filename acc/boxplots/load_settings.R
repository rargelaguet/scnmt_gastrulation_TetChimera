################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera/settings.R")
} else {
  stop("Computer not recognised")
}
io$diff.acc <- paste0(io$basedir,"/acc/results/differential/feature_level/lineages")
io$outdir <- paste0(io$basedir,"/acc/results/boxplots")

# sample_metadata[pass_accQC==T,.N,by="lineage"]

####################
## Define options ##
####################

# Define lineages
opts$lineages <- c(
  "Prelineage",
  # "Morula",
  "ICM",
  # "Epiblast",
  # "PrE",
  "TE_mural",
  "TE_polar"
)

# Define genomic contexts for accessibility
opts$acc.annos <- c(
  # "atac_peaks_2cell_2pn",
  # "atac_peaks_2cell_3pn",
  # "atac_peaks_4cell_3pn",
  # "atac_peaks_8cell_2pn",
  # "atac_peaks_8cell_3pn",
  # "atac_peaks_icm_2pn"
  # "ICM_H3K27ac",
  # "ICM_H3K4me3",
  # "TE_H3K27me3",
  # "distal_ICM_H3K27ac",
  # "distal_ICM_H3K4me3"
  "prom_200_200",
  # "genebody",
  # "prom_200_200_cgi",
  # "prom_200_200_noncgi"
  "LINE",
  "L1"
  # "CR1"
)
# opts$acc.annos <- NULL
# if (is.null(opts$acc.annos)) {
#   opts$acc.annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
# }

# Options for selecting differential hits
opts$fdr.threshold <- 0.10
opts$min.acc.diff <- 25


# Update metadata
sample_metadata <- sample_metadata %>% 
  .[lineage%in%opts$lineages] %>%
  .[pass_accQC==TRUE] %>%
  droplevels

table(sample_metadata$lineage)
