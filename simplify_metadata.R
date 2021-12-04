foo <- fread("/Users/argelagr/data/tet_chimera_nmtseq/backup/sample_metadata.txt")



foo[!is.na(cg_files),id_met:=cg_files %>% strsplit("/") %>% map_chr(9) %>% gsub(".tsv.gz","",.)]
foo[!is.na(gc_files),id_acc:=gc_files %>% strsplit("/") %>% map_chr(9) %>% gsub(".tsv.gz","",.)]
View(foo)

bar <- foo[,c("cell","id_rna","id_met","id_acc","plate","id_rna","method","embryo","ko_type","stage","class","tdTOM", "KDR-Cy7", "CD41-BV421")]
fwrite(bar,"/Users/argelagr/data/tet_chimera_nmtseq/sample_metadata.txt", quote=F, col.names = T, sep="\t", na="NA")



sce <- readRDS("/Users/argelagr/data/tet_chimera_nmtseq/processed/rna/SingleCellExperiment.rds")

metadata <- colData(sce) %>% as.data.table(keep.rownames = T) %>% setnames("rn","id_rna")

# where are plates1 and 2 for the E7.5???
plate2sample <- c(
  "E7.5_tet_chimera_plate3" = "E7.5_TET_TKO_KDR+",
  "E7.5_tet_crispr_plate5" = "E7.5_TET_TKO_crispr_rep1",
  "E7.5_tet_crispr_plate6" = "E7.5_TET_TKO_crispr_rep2",
  "E8.5_oct20_plate1" = "E8.5_WT_rep1",
  "E8.5_oct20_plate7" = "E8.5_WT_rep2",
  "E8.5_oct20_plate5" = "E8.5_WT_CD41+_rep1",
  "E8.5_oct20_plate3" = "E8.5_WT_CD41+_rep2",
  "E8.5_oct20_plate4" = "E8.5_TET_TKO_KDR+_CD41+_rep1",
  "E8.5_oct20_plate6" = "E8.5_TET_TKO_KDR+_CD41+_rep2",
  "E8.5_oct20_plate2" = "E8.5_TET_TKO_KDR+_rep1",
  "E8.5_oct20_plate8" = "E8.5_TET_TKO_KDR+_rep2"
)

metadata[,sample:=stringr::str_replace_all(plate,plate2sample)]
table(metadata$sample)
fwrite(metadata,"/Users/argelagr/data/tet_chimera_nmtseq/processed/sample_metadata_sce.txt", quote=F, col.names = T, sep="\t", na="NA")
