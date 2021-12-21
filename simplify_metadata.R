foo <- fread("/Users/argelagr/data/tet_chimera_nmtseq/sample_metadata.txt.gz")
foo <- fread("/Users/argelagr/data/tet_chimera_nmtseq/backup/sample_metadata.txt")



foo[!is.na(cg_files),id_met:=cg_files %>% strsplit("/") %>% map_chr(9) %>% gsub(".tsv.gz","",.)]
foo[!is.na(gc_files),id_acc:=gc_files %>% strsplit("/") %>% map_chr(9) %>% gsub(".tsv.gz","",.)]
View(foo)

bar <- foo[,c("cell","id_rna","id_met","id_acc","plate","id_rna","method","embryo","ko_type","stage","class")]
fwrite(bar,"/Users/argelagr/data/tet_chimera_nmtseq/sample_metadata.txt", quote=F, col.names = T, sep="\t", na="NA")



sce <- readRDS("/Users/argelagr/data/tet_chimera_nmtseq/processed/rna/SingleCellExperiment.rds")

metadata <- colData(sce) %>% as.data.table(keep.rownames = T) %>% setnames("rn","id_rna")

# create sample column
plate2sample <- c(
  "E7.5_tet_chimera_plate3" = "E7.5_TET_TKO",
  "E7.5_tet_crispr_plate5" = "E7.5_TET_TKO_crispr",
  "E7.5_tet_crispr_plate6" = "E7.5_TET_TKO_crispr",
  "E8.5_oct20_plate1" = "E8.5_WT_CD41+",
  "E8.5_oct20_plate2" = "E8.5_TET_TKO_CD41+",
  "E8.5_oct20_plate3" = "E8.5_WT_KDR+",
  "E8.5_oct20_plate4" = "E8.5_TET_TKO_KDR+",
  "E8.5_oct20_plate5" = "E8.5_WT_KDR+_CD41+",
  "E8.5_oct20_plate6" = "E8.5_TET_TKO_KDR+_CD41+",
  "E8.5_oct20_plate7" = "E8.5_WT",
  "E8.5_oct20_plate8" = "E8.5_TET_TKO"
)

metadata[,sample:=stringr::str_replace_all(plate,plate2sample)]
table(metadata$sample)

# create class column
# where are plates1 and 2 for the E7.5???
# plate2class <- c(
#   "E7.5_tet_chimera_plate3" = "E7.5_TET_TKO",
#   "E7.5_tet_crispr_plate5" = "E7.5_TET_TKO_crispr",
#   "E7.5_tet_crispr_plate6" = "E7.5_TET_TKO_crispr",
#   "E8.5_oct20_plate1" = "E8.5_WT_CD41+",
#   "E8.5_oct20_plate2" = "E8.5_TET_TKO_CD41+",
#   "E8.5_oct20_plate3" = "E8.5_WT_KDR+",
#   "E8.5_oct20_plate4" = "E8.5_TET_TKO_KDR+",
#   "E8.5_oct20_plate5" = "E8.5_WT_KDR+_CD41+",
#   "E8.5_oct20_plate6" = "E8.5_TET_TKO_KDR+_CD41+",
#   "E8.5_oct20_plate7" = "E8.5_WT",
#   "E8.5_oct20_plate8" = "E8.5_TET_TKO"
# )
# 
# metadata[,class:=stringr::str_replace_all(plate,plate2class)]
# table(metadata$class)

# create id_met and id_acc
metadata[!is.na(cg_files),id_met:=cg_files %>% gsub(".tsv.gz","",.)]
metadata[!is.na(gc_files),id_acc:=gc_files %>% gsub(".tsv.gz","",.)]

# remove unused columns
metadata[,c("cg_files","gc_files","sizeFactor","sierra_dna","sierra_rna","CD41.BV421","KDR.Cy7","tdTOM","i5","i7","well","column","class"):=NULL]

fwrite(metadata,"/Users/argelagr/data/tet_chimera_nmtseq/processed/sample_metadata_sce.txt.gz", quote=F, col.names = T, sep="\t", na="NA")



# foo <- fread("/Users/argelagr/data/tet_chimera_nmtseq/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# bar <- fread("/Users/argelagr/data/tet_chimera_nmtseq/processed/sample_metadata_sce.txt.gz") %>% .[,c("plate","sample")] %>% unique
# 
# baz <- foo %>% .[,sample:=NULL] %>% merge(bar,by="plate")
# table(baz$sample)
# fwrite(baz,"/Users/argelagr/data/tet_chimera_nmtseq/results/rna/mapping/sample_metadata_after_mapping.txt.gz", quote = F, col.names = T, na="NA", sep="\t")

