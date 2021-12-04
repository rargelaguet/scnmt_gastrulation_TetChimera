################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera//settings.R")
} else {
  stop("Computer not recognised")
}
io$diff.met <- paste0(io$basedir,"/met/results/differential/feature_level/lineages")
io$outdir <- paste0(io$basedir,"/met/results/boxplots")

####################
## Define options ##
####################

# Define lineages
opts$lineages <- c(
  "Prelineage",
  "Morula",
  "ICM",
  "Epiblast",
  "PrE",
  "TE_mural",
  "TE_polar"
)

# Define genomic contexts for methylation
opts$met.annos <- c(
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
  "prom_2000_2000"
  # "prom_2000_2000_cgi",
  # "prom_2000_2000_noncgi"
  # "LINE",
  # "L1"
  # "CR1"
)
# opts$met.annos <- NULL
# if (is.null(opts$met.annos)) {
#   opts$met.annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
# }

# Options for selecting differential hits
opts$fdr.threshold <- 0.10
opts$min.met.diff <- 25


# Update metadata
sample_metadata <- sample_metadata %>% 
  .[lineage%in%opts$lineages] %>%
  .[pass_metQC==TRUE]

table(sample_metadata$lineage)
