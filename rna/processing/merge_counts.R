library(data.table)
library(purrr)

matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

##############
## Metadata ##
##############

metadata <- fread("/Users/argelagr/data/tet_chimera_nmtseq/processed/sample_metadata.txt.gz")
metadata[,2] <- NULL
metadata %>% .[is.na(`CD41-BV421`), `CD41-BV421`:=FALSE] %>% .[is.na(`KDR-Cy7`), `KDR-Cy7`:=FALSE]

metadata_old_batch1 <- fread("/Users/argelagr/data/tet_chimera_nmtseq/processed/rna/tmp/batch1_metadata.txt") %>%
  .[,c("id","id_rna","plate","ko","stage")] %>%
  .[,c("KDR-Cy7","CD41-BV421","tdTOM"):=as.logical(FALSE)] %>%
  .[ko=="tet_tko", tdTOM:=TRUE] %>%
  .[ko=="host", tdTOM:=FALSE] %>% 
  .[,ko:=NULL] %>% .[,ko_type:="chimera"] %>%
  setnames("id","cell") %>%
  .[,class:=ifelse(tdTOM=="TRUE","E7.5_chimera_KO","E7.5_host_WT")] %>%
  .[,method:="rna"] %>%
  .[,embryo:=plate] %>%
  .[,c("cell","id_rna", "plate", "embryo", "stage", "method","ko_type", "class","tdTOM", "CD41-BV421", "KDR-Cy7")]

metadata_old_batch2 <- fread("/Users/argelagr/data/tet_chimera_nmtseq/processed/rna/tmp/batch2_metadata.gz") %>%
  .[,c("id_rna","plate","stage","tdTOM","CD41","KDR")] %>%
  .[,cell:=gsub("_GRCm38_hisat2_bam","",id_rna)] %>%
  setnames("KDR","KDR-Cy7") %>% setnames("CD41","CD41-BV421") %>%
  .[is.na(`CD41-BV421`), `CD41-BV421`:=FALSE] %>% .[is.na(`KDR-Cy7`), `KDR-Cy7`:=FALSE] %>%
  .[,ko_type:="chimera"] %>%
  .[,class:=ifelse(tdTOM=="TRUE","chimera_KO","host_WT")] %>%
  .[stage=="E7.5",class:=paste0("E7.5_",class)] %>% .[stage=="E8.5",class:=paste0("E8.5_",class)] %>%
  .[,method:="rna"] %>%
  .[,embryo:=plate] %>%
  .[,c("cell","id_rna", "plate", "embryo","stage", "method", "ko_type", "class","tdTOM", "CD41-BV421", "KDR-Cy7")]


metadata_old_batch2[plate=="tet_chimera_march20_plate4" & tdTOM=="TRUE",plate:="tet_chimera_march20_plate4_tet_tko"]
metadata_old_batch2[plate=="tet_chimera_march20_plate4" & tdTOM=="FALSE",plate:="tet_chimera_march20_plate4_wt"]
metadata_old_batch2[plate=="tet_chimera_march20_plate5" & tdTOM=="TRUE",plate:="tet_chimera_march20_plate5_tet_tko"]
metadata_old_batch2[plate=="tet_chimera_march20_plate5" & tdTOM=="FALSE",plate:="tet_chimera_march20_plate5_wt"]
metadata_old_batch2[plate=="tet_chimera_march20_plate6" & tdTOM=="TRUE",plate:="tet_chimera_march20_plate6_tet_tko"]
metadata_old_batch2[plate=="tet_chimera_march20_plate6" & tdTOM=="FALSE",plate:="tet_chimera_march20_plate6_wt"]
table(metadata_old_batch2$plate)

############
## Counts ##
############

counts_old_batch1.mtx <- fread("/Users/argelagr/data/tet_chimera_nmtseq/processed/rna/tmp/batch1_counts.tsv.gz") %>% matrix.please
counts_old_batch2.mtx <- fread("/Users/argelagr/data/tet_chimera_nmtseq/processed/rna/tmp/batch2_counts.gz") %>% matrix.please
counts.mtx <- fread("/Users/argelagr/data/tet_chimera_nmtseq/processed/rna/counts.tsv.gz") %>% matrix.please

# genes <- Reduce(intersect, list(rownames(counts.mtx),rownames(counts_old_batch1.mtx), rownames(counts_old_batch2.mtx)))
all(rownames(counts_old_batch1.mtx)==rownames(counts_old_batch2.mtx))
all(rownames(counts.mtx)==rownames(counts_old_batch2.mtx))

##################
## Merge counts ##
##################

counts_merged.mtx <- do.call("cbind",list(counts.mtx, counts_old_batch1.mtx, counts_old_batch2.mtx))

#####################
## Merge metadata ##
####################

metadata_old <- rbind(metadata_old_batch1,metadata_old_batch2) %>%
  .[,c("id_met","id_acc"):=as.character(NA)]

colnames(metadata_old)[!colnames(metadata_old)%in%colnames(metadata)]
colnames(metadata)[!colnames(metadata)%in%colnames(metadata_old)]
stopifnot(sort(colnames(metadata))==sort(colnames(metadata_old)))
# metadata_merged <- merge()

cols <- c("cell", "id_rna", "id_met", "id_acc", "method", "stage", "plate", "embryo", "ko_type", "class", "tdTOM", "KDR-Cy7", "CD41-BV421")
metadata_merged <- rbind(metadata[,..cols], metadata_old[,..cols])
stopifnot(colnames(counts_merged.mtx) %in% metadata_merged$id_rna)

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
  "E8.5_oct20_plate8" = "E8.5_TET_TKO",
  "tet_chimera_march20_plate1" = "E8.5_WT_CD41+",
  "tet_chimera_march20_plate2" = "E8.5_WT_KDR+",
  "tet_chimera_march20_plate3" = "E8.5_TET_TKO_KDR+",
  "tet_chimera_march20_plate4_tet_tko" = "E8.5_TET_TKO",
  "tet_chimera_march20_plate4_wt" = "E8.5_WT",
  "tet_chimera_march20_plate5_tet_tko" = "E7.5_TET_TKO",
  "tet_chimera_march20_plate5_wt" = "E7.5_WT",
  "tet_chimera_march20_plate6_tet_tko" = "E7.5_TET_TKO",
  "tet_chimera_march20_plate6_wt" = "E7.5_WT",
  "Tet-tko-chimera_embryo01" = "E7.5_TET_TKO",
  "Tet-tko-chimera_embryo02" = "E7.5_TET_TKO",
  "Tet-tko-chimera_embryo04" = "E7.5_WT",
  "Tet-tko-chimera_embryo05" = "E7.5_WT"
)
unique(metadata_merged$plate)[!unique(metadata_merged$plate)%in%names(plate2sample)]
names(plate2sample)[!names(plate2sample)%in%unique(metadata_merged$plate)]

stopifnot(sort(names(plate2sample))==sort(unique(metadata_merged$plate)))

metadata_merged[,sample:=stringr::str_replace_all(plate,plate2sample)]
table(metadata_merged$sample)
stopifnot(!is.na(metadata_merged$sample))

# TO-DO
metadata_merged[method=="mt",id_acc:=NA]

##########
## Save ##
##########

# Save plate metadata
plate_metadata <- metadata_merged[,.(N=.N),c("plate", "stage", "method", "ko_type", "class", "sample", "tdTOM", "KDR-Cy7", "CD41-BV421")]
fwrite(plate_metadata,"/Users/argelagr/data/tet_chimera_nmtseq/processed/plate_metadata.txt.gz", quote=F, col.names = T, sep="\t", na="NA")

# Save cell metadata
fwrite(metadata_merged,"/Users/argelagr/data/tet_chimera_nmtseq/processed/sample_metadata_merged.txt.gz", quote=F, col.names = T, sep="\t", na="NA")

write.table(counts_merged.mtx, "/Users/argelagr/data/tet_chimera_nmtseq/processed/rna_new/counts.txt", sep="\t", col.names=T, row.names=T, quote=F, na="NA")

##########
## TEST ##
##########

# foo <- fread("/Users/argelagr/data/not_used/scnmt_gastrulation_TetChimera/metadata.txt.gz") %>%
#   .[,c("id_rna","plate","stage","tdTOM","CD41","KDR")]
# all(foo$id_rna == metadata_old_batch2$id_rna)
# 
# bar <- fread("/Users/argelagr/data/not_used/scnmt_gastrulation_TetChimera/TO_MERGE_scnmt_gastrulation_TetKO/metadata.txt.gz")
# all(bar$id_rna == metadata_old_batch1$id_rna)
