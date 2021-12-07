foo <- fread("/Users/argelagr/data/gastrulation_multiome_10x/processed_new/atac/archR/PeakCalls/peaks_archR_macs2.bed.gz") %>%
  setnames(c("chr","start","end")) %>%
  .[,id:=sprintf("%s:%s-%s",chr,start,end)] %>% 
  .[,strand:="*"] %>%
  .[,c("chr","start","end","strand","id")]


fwrite(foo, "/Users/argelagr/data/tet_chimera_nmtseq/features/genomic_contexts/multiome_peaks.bed.gz", col.names = F, sep="\t", quote=F)
