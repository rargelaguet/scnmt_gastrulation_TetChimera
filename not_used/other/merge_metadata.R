library(data.table)
library(purrr)

########
## WT ##
########

cols <- c("sample", "id_rna", "id_met", "id_acc", "pass_rnaQC", "pass_metQC", "pass_accQC", "stage", "lineage10x","embryo","plate")
a <- fread("/Users/ricard/data/scnmt_gastrulation_merged_dnmtKO/sample_metadata_wt.txt") %>% .[,..cols]
a[,genotype:="WT"]
a[,lineage10x:=stringr::str_replace_all(lineage10x,"_"," ")]

########
## KO ##
########

cols <- c("sample", "id_rna", "id_met", "id_acc", "pass_rnaQC", "pass_metQC", "pass_accQC", "stage", "lineage10x","embryo")
b <- fread("/Users/ricard/data/scnmt_gastrulation_merged_dnmtKO/sample_metadata_ko.txt") %>% 
  .[,..cols]
b[,plate:=embryo]
b[,genotype:="KO"]

###########
## Merge ##
###########

ab <- plyr::rbind.fill(a,b) %>% as.data.table

fwrite(ab, "/Users/ricard/data/scnmt_gastrulation_merged_dnmtKO/sample_metadata.txt", sep="\t", col.names=T, row.names=F, na="NA", quote=F)

# write.table(a_filt, file="/Users/ricard/data/gastrulation_merged2/rna/E3.5_counts.txt", sep="\t", col.names=T, row.names=T, quote=F)
