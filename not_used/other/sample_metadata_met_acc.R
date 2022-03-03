foo <- fread("/Users/argelagr/data/tet_chimera_nmtseq/backup/sample_metadata.txt")
foo[,c("cell","cg_files","gc_files")] %>% .[!is.na(cg_files)]

all(sample_metadata$cell%in%foo$cell)

metadata[!is.na(cg_files),id_met:=cg_files %>% gsub(".tsv.gz","",.)]
metadata[!is.na(gc_files),id_acc:=gc_files %>% gsub(".tsv.gz","",.)]

foo[!is.na(cg_files),id_met:=cg_files %>% strsplit("/") %>% map_chr(9) %>% gsub(".tsv.gz","",.)]
foo[!is.na(gc_files),id_acc:=gc_files %>% strsplit("/") %>% map_chr(9) %>% gsub(".tsv.gz","",.)]


files <- list.files("/Users/argelagr/data/tet_chimera_nmtseq/processed/met/cpg_level", pattern = "*.tsv.gz") %>% gsub(".tsv.gz","",.)
mean(files %in% foo$id_met)

files <- list.files("/Users/argelagr/data/tet_chimera_nmtseq/processed/acc/gpc_level", pattern = "*.tsv.gz") %>% gsub(".tsv.gz","",.)
mean(files %in% foo$id_acc)

bar <- foo[,c("cell","id_met","id_acc")]
fwrite(bar,"/Users/argelagr/data/tet_chimera_nmtseq/sample_metadata_metacc.txt.gz", quote=F, col.names = T, sep="\t", na="NA")






foo <- fread("/Users/argelagr/data/tet_chimera_nmtseq/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
bar <- fread("/Users/argelagr/data/tet_chimera_nmtseq/sample_metadata_metacc.txt.gz")
all(foo$cell%in%bar$cell)

baz <- foo %>% .[,c("id_met","id_acc"):=NULL] %>% merge(bar,by="cell")
mean(!is.na(baz$id_met))
mean(!is.na(baz$id_acc))
fwrite(baz,"/Users/argelagr/data/tet_chimera_nmtseq/results/rna/mapping/sample_metadata_after_mapping.txt.gz", quote = F, col.names = T, na="NA", sep="\t")
