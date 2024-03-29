library(data.table)
library(purrr)

sample_metadata <- fread("/bi/group/reik/ricard/data/tet_chimera_nmtseq/results_new/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")

samples.to.remove <- sample_metadata[method=="mt" & !is.na(id_met),id_met]
for (i in samples.to.remove) {
	files.to.remove <- grep(i,list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/met/feature_level/tmp",full.names=T),value=T)
	file.remove(files.to.remove)
}


samples.to.remove <- sample_metadata[method=="mt" & !is.na(id_met),id_met] %>% gsub("_CpG","",.)
for (i in samples.to.remove) {
	files.to.remove <- grep(i,list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/acc/feature_level/tmp",full.names=T),value=T)
	file.remove(files.to.remove)
}

