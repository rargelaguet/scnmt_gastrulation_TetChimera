suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/scnmt_gastrulation_TetChimera"
  io$atlas.basedir <- "/Users/ricard/data/gastrulation10x"
  # io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else if (grepl("rargelaguet.local",Sys.info()['nodename'])) {
  io$basedir <- "/Users/rargelaguet/data/scnmt_gastrulation_TetChimera"
  io$atlas.basedir <- "/Users/rargelaguet/data/gastrulation10x"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_gastrulation_TetChimera"
  io$atlas.basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
  # io$gene_metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else if (Sys.info()[['nodename']]=="BI2404M") {
  io$basedir <- "/Users/argelagr/data/tet_chimera_nmtseq"
  io$atlas.basedir <- "/Users/argelagr/data/gastrulation10x"
  io$atlas.basedir <- "/Users/argelagr/data/gastrulation10x"
  io$multiome.basedir <- "/Users/argelagr/data/gastrulation_multiome_10x"
} else if (grepl("bi2228m", tolower(Sys.info()["nodename"]))){
  io$basedir <- "/Users/clarks/data/nmt_tet_chimera"
} else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
  if (grepl("Clark", Sys.info()['effective_user'])) {
    io$basedir       <- "/bi/scratch/Stephen_Clark/multiome/resilio"
  } else if (grepl("argelag", Sys.info()['effective_user'])) {
    io$basedir <- "/bi/group/reik/ricard/data/tet_chimera_nmtseq"
    io$atlas.basedir <- "/bi/group/reik/ricard/data/pijuansala2019_gastrulation10x"
  	io$multiome.basedir <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x"
  }
} else {
  stop("Computer not recognised")
}

# Metadata
io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
io$gene_metadata <- paste0(io$basedir, "/features/gene_metadata/Mmusculus_genes_BioMart.87.txt")
io$features.dir <- paste0(io$basedir, "/features/genomic_contexts")

# RNA
# io$rna.sce <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")
io$rna.sce <- paste0(io$basedir,"/processed/rna/SingleCellExperiment.rds")
# io$rna.counts <- paste0(io$basedir,"/rna/counts.tsv.gz")
# io$rna.diff <- paste0(io$basedir,"/rna/results/differential/lineages")
# io$rna.stats <- paste0(io$basedir,"/rna/results/stats/rna_stats.txt.gz")


# Methylation
io$met_data_raw <- paste0(io$basedir,"/processed/met/cpg_level")
io$met_data_raw_pseudobulk <- paste0(io$met_data_raw,"/pseudobulk")
io$met_data_parsed <- paste0(io$basedir,"/processed/met/feature_level")
io$met_data_parsed_pseudobulk <- paste0(io$met_data_parsed,"/pseudobulk")
io$met.stats <- paste0(io$basedir,"/met/results/stats/sample_stats.txt")
io$met.diff <- paste0(io$basedir,"/met/results/differential/feature_level/lineages")

# Accessibility
io$acc_data_raw <- paste0(io$basedir,"/processed/acc/gpc_level")
io$acc_data_parsed <- paste0(io$basedir,"/processed/acc/feature_level")
io$acc_data_raw_pseudobulk <- paste0(io$acc_data_raw,"/pseudobulk")
io$acc_data_parsed_pseudobulk <- paste0(io$acc_data_parsed,"/pseudobulk")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/sample_stats.txt")
io$acc.diff <- paste0(io$basedir,"/acc/results/differential/feature_level/lineages")

# 10x multiome information
io$multiome.metadata <- paste0(io$multiome.basedir,"/sample_metadata.txt.gz")
io$multiome.differential <- paste0(io$multiome.basedir,"/results/differential")
io$multiome.sce <- paste0(io$multiome.basedir,"/processed/SingleCellExperiment.rds")
io$multiome.marker_peaks <- paste0(io$multiome.basedir,"/results/atac/archR/differential/PeakMatrix/markers/marker_peaks.txt.gz")  # this needs to be updated

# PijuanSala2019 information
io$atlas.metadata <- paste0(io$atlas.basedir,"/sample_metadata.txt.gz")
io$atlas.marker_genes <- paste0(io$atlas.basedir,"/results/marker_genes/all_stages/marker_genes.txt.gz")
io$atlas.differential <- paste0(io$atlas.basedir,"/results/differential")
io$atlas.average_expression_per_celltype <- paste0(io$atlas.basedir,"/results/marker_genes/avg_expr_per_celltype_and_gene.txt.gz")
io$atlas.sce <- paste0(io$atlas.basedir,"/processed/SingleCellExperiment.rds")


#############
## Options ##
#############

opts <- list()

opts$celltypes = c(
	"Epiblast",
	"Primitive_Streak",
	"Caudal_epiblast",
	"PGC",
	"Anterior_Primitive_Streak",
	"Notochord",
	"Def._endoderm",
	"Gut",
	"Nascent_mesoderm",
	"Mixed_mesoderm",
	"Intermediate_mesoderm",
	"Caudal_Mesoderm",
	"Paraxial_mesoderm",
	"Somitic_mesoderm",
	"Pharyngeal_mesoderm",
	"Cardiomyocytes",
	"Allantois",
	"ExE_mesoderm",
	"Mesenchyme",
	"Haematoendothelial_progenitors",
	"Endothelium",
	"Blood_progenitors",
	# "Blood_progenitors_2",
	# "Blood_progenitors_2",
	# "Erythroid1",
	# "Erythroid2",
	# "Erythroid3",
	"early_Erythroid",
	"late_Erythroid",
	"NMP",
	"Neurectoderm",
	# "Rostral_neurectoderm",
	# "Caudal_neurectoderm",
	"Neural_crest",
	"Forebrain_Midbrain_Hindbrain",
	"Spinal_cord",
	"Surface_ectoderm",
	"Visceral_endoderm",
	"ExE_endoderm",
	"ExE_ectoderm",
	"Parietal_endoderm"
)

opts$celltype.colors = c(
	"Epiblast" = "#635547",
	"Primitive_Streak" = "#DABE99",
	"Caudal_epiblast" = "#9e6762",
	"PGC" = "#FACB12",
	"Anterior_Primitive_Streak" = "#c19f70",
	"Notochord" = "#0F4A9C",
	"Def._endoderm" = "#F397C0",
	"Gut" = "#EF5A9D",
	"Nascent_mesoderm" = "#C594BF",
	"Mixed_mesoderm" = "#DFCDE4",
	"Intermediate_mesoderm" = "#139992",
	"Caudal_Mesoderm" = "#3F84AA",
	"Paraxial_mesoderm" = "#8DB5CE",
	"Somitic_mesoderm" = "#005579",
	"Pharyngeal_mesoderm" = "#C9EBFB",
	"Cardiomyocytes" = "#B51D8D",
	"Allantois" = "#532C8A",
	"ExE_mesoderm" = "#8870ad",
	"Mesenchyme" = "#cc7818",
	"Haematoendothelial_progenitors" = "#FBBE92",
	"Endothelium" = "#ff891c",
	"Blood_progenitors_1" = "#f9decf",
	"Blood_progenitors_2" = "#c9a997",
	"Blood_progenitors" = "#c9a997",
	"Erythroid1" = "#C72228",
	"Erythroid2" = "#f79083",
	"Erythroid3" = "#EF4E22",
	"late_Erythroid" = "#C72228",
	"early_Erythroid" = "#EF4E22",
	"NMP" = "#8EC792",
	"Rostral_neurectoderm" = "#65A83E",
	"Neurectoderm" = "#65A83E",
	"Caudal_neurectoderm" = "#354E23",
	"Neural_crest" = "#C3C388",
	"Forebrain_Midbrain_Hindbrain" = "#647a4f",
	"Spinal_cord" = "#CDE088",
	"Surface_ectoderm" = "#f7f79e",
	"Visceral_endoderm" = "#F6BFCB",
	"ExE_endoderm" = "#7F6874",
	"ExE_ectoderm" = "#989898",
	"Parietal_endoderm" = "#1A1A1A"
)

opts$plates <- c(
  "E7.5_tet_chimera_plate3",
  "E7.5_tet_crispr_plate5",
  "E7.5_tet_crispr_plate6",
  "E8.5_oct20_plate1",
  "E8.5_oct20_plate2",
  "E8.5_oct20_plate3",
  "E8.5_oct20_plate4",
  "E8.5_oct20_plate5",
  "E8.5_oct20_plate6",
  "E8.5_oct20_plate7",
  "E8.5_oct20_plate8",
  "Tet-tko-chimera_embryo01",
  "Tet-tko-chimera_embryo02",
  "Tet-tko-chimera_embryo04",
  "Tet-tko-chimera_embryo05",
  "tet_chimera_march20_plate1",
  "tet_chimera_march20_plate2",
  "tet_chimera_march20_plate3",
  "tet_chimera_march20_plate4_tet_tko",
  "tet_chimera_march20_plate4_wt",
  "tet_chimera_march20_plate5_tet_tko",
  "tet_chimera_march20_plate5_wt",
  "tet_chimera_march20_plate6_tet_tko",
  "tet_chimera_march20_plate6_wt"
)

opts$samples <- c(
  "E7.5_WT",
  "E7.5_TET_TKO",
  "E7.5_TET_TKO_crispr",
  "E8.5_WT",
  "E8.5_WT_CD41+",
  "E8.5_WT_KDR+",
  "E8.5_WT_KDR+_CD41+",
  "E8.5_TET_TKO",
  "E8.5_TET_TKO_KDR+",
  "E8.5_TET_TKO_CD41+",
  "E8.5_TET_TKO_KDR+_CD41+"
)

opts$chr <- paste0("chr",c(1:19,"X","Y"))

opts$context.colors <- c("CG"="#F8766D", "GC"="#00BFC4")

opts$class.colors <- c("WT" = "#4F94CD", "Tet-TKO" = "#EE4000", "TET-TKO" = "#EE4000")
