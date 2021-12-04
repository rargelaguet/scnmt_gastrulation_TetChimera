
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera/settings.R")
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/acc/dimensionality_reduction")

####################
## Define options ##
####################

# Define which annotations to look at
opts$annos <- c(
  # "prom_200_200" = "Promoters"
  "H1_distal_H3K27ac" = "hESCs distal H3K27ac",
  "H1_H3K4me3" = "hESCs H3K4me3",
  "H1_H3K4me1" = "hESCs H3K4me1",
  # "CGI" = "CpG islands"
  # "trophoblast_cell_ChIP-seq_H3K27me3-human" = "Trophoblast H3K27me3",
  # "trophoblast_cell_ChIP-seq_H3K4me1-human" = "Trophoblast H3K4me1"
)

opts$stages <- c(
  # "E3", 
  # "E4", 
  "E5",
  "E6"
  # "E7",
  # "E10"
)

opts$labs <- c(
  # "lanner",
  "nichols"
)


# Filtering options
opts$min.GpCs <- 5          # minimum number of CpG sites per feature and cell
opts$min.coverage <- 0.20   # minimum coverage (fraction of cells with at least min.GpC measurements)
opts$nfeatures <- 1500     # number of features per view (filter based on variance)

# Output file
io$outfile = sprintf("%s/hdf5/model_%s_%s.hdf5",io$outdir,paste(names(opts$annos), collapse="_"), paste(opts$stages, collapse="_"))

############################
## Update sample metadata ##
############################

sample_metadata <- sample_metadata %>% 
  .[pass_accQC==T & stage%in%opts$stages & lab%in%opts$labs] %>%
  # .[,c("id_acc","id_rna","stage","lab")] %>%
  merge(fread(io$acc.stats)[,c("id_acc","mean","coverage")], by="id_acc")

