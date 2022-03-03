library(data.table)
library(purrr)

##################
## scNMT-seq met ##
##################

indir <- "/Users/argelagr/data/scnmt_gastrulation/met/cpg_level"
outdir <- "/Users/argelagr/data/scnmt_gastrulation_argelaguet2019/processed/met/cpg_level/new"

files <- list.files(indir, full.names=F, pattern="*.gz")

for (i in files) {
	print(sprintf("%s (%s/%s)",i,match(i,files),length(files)))
	tmp <- fread(file.path(indir,i), select=c(1,2,5)) %>% setnames(c("chr","pos","rate"))
	fwrite(tmp, file.path(outdir,i), sep="\t", col.names=T, quote=F, na="NA")
}


###################
## scNMT-seq acc ##
##################3

indir <- "/Users/argelagr/data/scnmt_gastrulation/acc/gpc_level"
outdir <- "/Users/argelagr/data/scnmt_gastrulation_argelaguet2019/processed/acc/gpc_level"

files <- list.files(indir, full.names=F, pattern="*.gz")

for (i in files) {
	print(sprintf("%s (%s/%s)",i,match(i,files),length(files)))
	tmp <- fread(file.path(indir,i), select=c(1,2,5)) %>% setnames(c("chr","pos","rate"))
	fwrite(tmp, file.path(outdir,i), sep="\t", col.names=T, quote=F, na="NA")
}