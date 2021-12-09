###########################
## Load default settings ##
###########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/human_embryo_multiomics/settings.R")
  source("/Users/ricard/human_embryo_multiomics/met/differential/feature_level/lineages/analysis/utils.R")
} else {
  stop("Computer not recognised")
}

io$indir <- paste0(io$basedir,"/met/results/differential/feature_level/lineages")
io$outdir <- paste0(io$basedir,"/met/results/differential/feature_level/lineages/pdf")

# sample_metadata[lineage%in%c("TE_mural","TE_polar") & pass_metQC==TRUE,.N,by="lineage"]

#############
## Options ##
#############

opts$comparisons <- c(
  "Prelineage_vs_Morula",
  "Prelineage_vs_ICM",
  "Prelineage_vs_Epiblast",
  "Prelineage_vs_TE_mural",
  "Prelineage_vs_TE_polar",
  "Prelineage_vs_PrE",
  "Morula_vs_ICM",
  "Morula_vs_Epiblast",
  "Morula_vs_TE_mural",
  "Morula_vs_TE_polar",
  "Morula_vs_PrE",
  "ICM_vs_Epiblast",
  "ICM_vs_PrE",
  "ICM_vs_TE_mural",
  "ICM_vs_TE_polar",
  "Epiblast_vs_TE_mural",
  "Epiblast_vs_TE_polar",
  "Epiblast_vs_PrE",
  "PrE_vs_TE_mural",
  "PrE_vs_TE_polar",
  "TE_mural_vs_TE_polar"
)

# Select genomic contexts
# opts$annos <- NULL
# if (is.null(opts$annos)) {
#   opts$annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
# }
opts$annos <- c(
  # "2cell_H3K27me3",
  # "4cell_H3K27me3",
  # "8cell_H3K27ac",
  "Alu",
  # "atac_peaks_2cell_2pn",
  # "atac_peaks_4cell_3pn",
  # "atac_peaks_8cell_2pn",
  # "atac_peaks_icm_2pn",
  "CR1",
  "CGI",
  # "distal_4cell_3PN_H3K4me3",
  # "distal_8cell_H3K27ac",
  # "distal_ICM_H3K27ac",
  # "distal_ICM_H3K4me3",
  "ERV1",
  "ERVL",
  "genebody",
  # "ICM_H3K27me3",
  # "ICM_H3K4me3",
  "L1",
  "L2",
  "LINE",
  # "prom_200_200_cgi",
  # "prom_200_200_noncgi",
  "prom_2000_2000_cgi",
  "prom_2000_2000_noncgi"
  # "TE_H3K27me3"
)

###############
## Load data ##
###############

# Load precomputed differential results
diff.results <- lapply(opts$comparisons, function(i) 
  lapply(opts$annos, function(j) {
    file <- sprintf("%s/%s_%s.txt.gz",io$indir,i,j)
    if (file.exists(file)) {
      fread(file) %>% .[,anno:=as.factor(j)]
    } else {
      cat(sprintf("%s does not exist\n",file))
    }
  }
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist


####################
## Define options ##
####################

# Minimum differential levels (%) for statistical significance
opts$min.diff <- 25

# Multiple testing correction
opts$threshold_fdr <- 0.10

#################
## Filter data ##
#################

# Define statistically significant hits
diff.results[,sig:=(padj_fdr<=opts$threshold_fdr & abs(diff)>opts$min.diff)]

# Filter genomic contexs with very few hits
diff.results <- diff.results[,N:=.N,by=c("anno","comparison")] %>% .[N>1000] %>% .[,N:=NULL]

