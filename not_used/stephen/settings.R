suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$atlas.basedir <- "/Users/ricard/data/gastrulation10x"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  
} else if (grepl("bi2228m", tolower(Sys.info()["nodename"]))){
  io$basedir <- "/Users/clarks/data/nmt_tet_chimera"
} else if (grepl("pebble", Sys.info()["nodename"])){
  io$basedir <- "/bi/scratch/Stephen_Clark/tet_chimera_nmtseq/"
  } else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt")

io$met_data_raw <- paste0(io$basedir,"/met/cpg_level")
io$met_data_raw.pseudobulk <- paste0(io$met_data_raw,"/pseudobulk")
io$met_data_parsed <- paste0(io$basedir,"/met/feature_level")
io$met_data_parsed.pseudobulk <- paste0(io$met_data_parsed,"/pseudobulk")
io$met.stats <- paste0(io$basedir,"/met/results/stats/sample_stats.txt")
io$met.diff <- paste0(io$basedir,"/met/results/differential/feature_level/lineages")


io$acc_data_raw <- paste0(io$basedir,"/acc/gpc_level")
io$acc_data_parsed <- paste0(io$basedir,"/acc/feature_level")
io$acc_data_raw.pseudobulk <- paste0(io$acc_data_raw,"/pseudobulk")
io$acc_data_parsed.pseudobulk <- paste0(io$acc_data_parsed,"/pseudobulk")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")
io$acc.diff <- paste0(io$basedir,"/acc/results/differential/feature_level/lineages")


# io$rna.sce <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")
io$rna.sce <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")
# io$rna.counts <- paste0(io$basedir,"/rna/counts.tsv.gz")
io$rna.diff <- paste0(io$basedir,"/rna/results/differential/lineages")
io$rna.stats <- paste0(io$basedir,"/rna/results/stats/rna_stats.txt.gz")
io$rna.pseudotime <- paste0(io$basedir,"/rna/results/pseudotime/pseudotime.tsv.gz")

io$features.dir <- paste0(io$basedir,"/features/filt")
io$features.tsv <- paste0(io$features.dir, "/multiome_peaks.bed.gz")
io$markerPeaks <- "/bi/scratch/Stephen_Clark/tet_chimera_nmtseq/features/filt/marker_peaks.txt.gz"
# io$cpg.density <- paste0(io$basedir,"/met/stats/features/cpg_density_perfeature.txt.gz")

#"/bi/scratch/Stephen_Clark/tet_chimera_nmtseq/features/gene_metadata/Mmusculus_genes_BioMart.87.txt"
io$gene_metadata <- paste0(io$basedir, "/features/gene_metadata/Mmusculus_genes_BioMart.87.txt")

# Atlas information
io$atlas.metadata <- paste0(io$atlas.basedir,"/sample_metadata.txt.gz")
io$atlas.marker_genes <- paste0(io$atlas.basedir,"/results/marker_genes/marker_genes.txt.gz")
io$atlas.differential <- paste0(io$atlas.basedir,"/results/differential")
io$atlas.average_expression_per_celltype <- paste0(io$atlas.basedir,"/results/marker_genes/avg_expr_per_celltype_and_gene.txt.gz")
io$atlas.sce <- paste0(io$atlas.basedir,"/processed/SingleCellExperiment.rds")




#############
## Options ##
#############

opts <- list()

# opts$stages <- c(
#   "E3", 
#   "E4", 
#   "E5",
#   "E6",
#   "E7"
#   # "E10"
# )
# 
# opts$lineages <- c(
#   "Prelineage",
#   "Morula",
#   "ICM",
#   "Epiblast",
#   "PrE",
#   "TE_mural",
#   "TE_polar"
# )
# 
# 
# opts$labs <- c(
#   "lanner", 
#   "nichols",
#   "Petropolous"
# )
# 
# opts$chr <- paste0("chr",c(1:22,"X","Y"))
# 
# opts$marker_genes <- list(
#   "ZGA" = c("DUXA", "ZSCAN5B", "ZSCAN4", "MBD3L2", "SCML2","TRIM43", "RFPL1", "PRAME"), # DUX4 (not expressed)
#   "Prelineage (early)" = c("KPNA7", "ACTL8", "NLRP13", "WEE2", "PCP4L1", 
#                            "NLRP4", "TCL1A", "KHDC3L", "OOEP", "PABPN1L", "H1FOO", "FAM46C", 
#                             "DCDC2", "GDF9", "NLRP11", "OOSP2", "PARL", "OTX2", 
#                            "ZNF280B", "RARRES2", "BMP15", "TLE6", "ZP2", "UNC13C", "GOLM1", 
#                            "RND1", "LRRTM4", "NLRP9", "HERC6", "ZP3", "PADI6", "PLEK2", 
#                            "NLRP5", "CCDC25", "HDC", "ACCSL", "HSPB1", "HABP2", "SERPINF1", 
#                            "TUBB8", "USP2", "PDE8B", "ARMC2", "YPEL5", "GDAP1", "C6orf52", 
#                            "LEUTX", "CNNM2", "HSF2BP", "ZP4", "CXorf67", "RFPL2", "BCL2L10", 
#                            "OSBPL10", "TBPL1", "PTTG1", "PRAMEF11", "MVP", "ZHX3", "TPRX1", 
#                            "HHLA2", "SPDL1", "MBD3L3", "ZSCAN5B", "CUL5", "DCAF4L1", "ZSCAN4", 
#                            "FAM13A", "KLF17", "SYTL3", "PRAMEF2", "CCNA1", "B3GNT2"), # DNMT1, UBXN11, DPPA5
#   "Prelineage (early&middle)" = c("DUXA", "TCL1A", "MFSD2A", "ZNF280A", "ULK2", "ACTL8", "ZSCAN5B", "RGS2", "SOX15", 
#                                   "SERTAD3", "KLF17", "TBPL1", "FOXR1", "ZIM3", "TMEM92", "DPPA3", "PTTG1", "PFKFB3", 
#                                   "YPEL5", "DPRX", "BIK", "DPPA5", "OOEP", "TPRX1", "ARGFX"),
#   "Prelineage (middle)" = c("DPRX", "BIK", "RNF11", "GABARAPL1", "ARGFX", "USP7", "RBP7", "PFKFB3"),
#   "Morula" = c("PPP1R14A", "PRAP1", "DHRS3", "CYP26A1", "ALPL", "ITM2C", "TMEM37"),
#   "ICM" = c("IFI16", "IL6R", "RAB15", "ZC3HAV1", "UACA", "ASRGL1", "PYGB", "PHLDA1"),   # NOT VERY WELL DEFINED.......
#   "ICM_Epiblast" = c("IFI16", "IL6R", "PIM2", "GDF3", "UACA", "ZC3HAV1", "GPX2", "RAB15", 
#                      "ETV4", "ASRGL1", "ALDH2", "IFITM1", "PHLDA1", "PYGB", "SPARC", "SAT1", "CBFA2T2"),
#   "Epiblast" = c("NODAL", "GDF3", "TDGF1", "IFITM1", "IFITM3", "ETV4", "PIM2", "LEFTY1", "LEFTY2", "FBP1"),
#   "PrE" = c("COL4A1","RSPO3","GATA4","SOX17","HNF1B","APOA1","P4HA1","TRIM2","IGF1","MYL4","COL4A2","P4HA1","FN1","SOX17","PDGFRA"),
#   # "TE (polar)" = c("GCM1", "GREM2", "SERPINB9", "PGF", "NRP1", "CSF3R", "SLCO4A1", "CCR7", "SP6", "CYP19A1", "PAPOLA", "S1PR2", "DLX5", "MUC15"),
#   "TE (polar)" = c("GCM1", "PGF", "CGA", "S1PR2", "CCR7", "GREM2", "LGALS3", "CBLB", 
#                    "MUC15", "SP6", "SDC1", "KRT23", "CCKBR", "CYP19A1", 
#                    "RHOBTB1", "PTN", "CYP11A1", "VGLL1", "NR2F2", "SLC38A1"),  # SERPINB9, DLG5
#   "TE (mural)" = c("S100A6", "FABP3", "LRP2", "GPRC5A", "ATP6V1B1"),
#   # "TE" = c("KRT19", "DLX3", "CCKBR", "PTN", "SLC38A1", "ABCG2", "CD24", "CYP11A1", "GATA2", "S100A6", "NR2F2", "VGLL1", "HAPLN1", "PWWP2B", "CA12", "DLG5")
#   "TE" = c("ABCG2", "KRT19", "DLX3", "GATA2", "SLC38A1", "STS", "PTN", 
#            "CA12", "CD24", "CCKBR", "S100A6", "CYP11A1", "FABP3", "RAB31", 
#            "LRP2", "HAPLN1", "MPP1", "SP6", "CLDN4", "VGLL1", "SMAGP", "NR2F2", 
#            "TACSTD2", "PRSS8", "RAB11FIP1", "EMP2", "GATA3", "DLG5"),
#   "TETs" = c("TET1","TET2","TET3"),
#   "DNMTs" = c("DNMT1","DNMT3A","DNMT3B","DNMT3L"),
#   "DPPAs" = c("DPPA2","DPPA3","DPPA4","DPPA5")
# )
# 
# ##########################
# ## Load sample metadata ##
# ##########################


metacols <- c(
  "cell",
  "stage_lineage", # needs to be defined
  "celltype.mapped",
  "class",
  "id_rna",
  "cg_files",
  "gc_files",
  "id_met", # needs to be defined
  "id_acc", # needs to be defined
  "tdTOM",
  "KDR-Cy7",
  "CD41-BV421",
  "markers",
  "pass_rnaQC",
  "pass_metQC",
  "pass_accQC"
)

sample_metadata <- fread(io$metadata)

# merge subclasses of celltypes
sample_metadata <- sample_metadata[, celltype.mapped := gsub("[1-9]", "", celltype.mapped)]%>% 
  .[,stage_lineage := paste(stage, celltype.mapped, sep = "_")]  
  

# add in met/acc id
sample_metadata[, id_met := gsub(".tsv.gz", "", basename(cg_files))]
sample_metadata[, id_acc := gsub(".tsv.gz", "", basename(gc_files))]


# condense flow sort info
sample_metadata[`KDR-Cy7` == TRUE & `CD41-BV421` == TRUE, markers := "double_positive"]
sample_metadata[`KDR-Cy7` == TRUE & `CD41-BV421` == FALSE, markers := "KDR_positive"]
sample_metadata[`KDR-Cy7` == FALSE & `CD41-BV421` == TRUE, markers := "CD41_positive"]


# subset to selected columns
sample_metadata <- sample_metadata[, ..metacols]

# 
# # factor.cols <- c("id_rna","id_met","id_acc","stage","lineage","lab","plate","embryo")
# # factor.cols <- c("id_rna","id_met","id_acc","stage","lab","plate","embryo")
# 
# sample_metadata <- fread(io$metadata) %>% 
#   .[stage%in%opts$stages & lineage%in%opts$lineages] %>%
#   .[,lineage:=factor(lineage,levels=opts$lineages)] %>%
#   .[,stage:=factor(stage,levels=opts$stages)]
#   # .[,(factor.cols):=lapply(.SD, as.factor),.SDcols=(factor.cols)] %>%
#   # droplevels
# 
# 
# 
# # sample_metadata[,lineage2:=lineage]
# # sample_metadata[lineage2%in%c("TE_mural","TE_polar"),lineage2:="TE"]
# # sample_metadata[lineage2%in%c("Morula"),lineage2:="Prelineage"]
# # sample_metadata[lineage2%in%c("Epiblast","ICM"),lineage2:="Epiblast_ICM"]
# # table(sample_metadata$lineage2)
# # fwrite(sample_metadata, io$metadata, sep="\t", na="NA", quote=F)
# 
# 
# 
# #####################
# ## Load file paths ##
# #####################
# 
# metadata <- fread(io$metadata)
# 
# cg_cells <- metadata[pass_metQC == TRUE, id_met]
# gc_cells <- metadata[pass_accQC == TRUE, id_acc]
# 
# cg_files <- dir(
#   io$met_data_raw,
#   pattern = ".tsv.gz$",
#   full = TRUE,
#   recursive = TRUE
# ) %>% 
#   set_names(basename(.) %>% gsub(".tsv.gz", "", .))
# 
# io$met_files_raw <- map(cg_cells, ~cg_files[names(cg_files) == .x])
# 
# gc_files <- dir(
#   io$acc_data_raw,
#   pattern = ".tsv.gz$",
#   full = TRUE,
#   recursive = TRUE
# ) %>% 
#   set_names(basename(.) %>% gsub(".tsv.gz", "", .))
# 
# io$acc_files_raw <- map(gc_cells, ~gc_files[names(gc_files) == .x])
# 
# 
# 
# 
# if (any(map_int(io$met_files_raw, length) != 1)){
#   warning("non unique or missing CpG files")
# }
# if (any(map_int(io$acc_files_raw, length) != 1)){
#   warning("non unique or missing GpC files")
# }

