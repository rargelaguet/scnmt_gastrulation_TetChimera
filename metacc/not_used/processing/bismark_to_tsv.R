library(data.table)
library(purrr)

##################
## scMT-seq met ##
##################

outdir <- "/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/met/cpg_level/new"
sample_metadata <- fread("/bi/group/reik/ricard/data/tet_chimera_nmtseq/results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")

files <- list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/original", full.names=T)

# Sanity check
files.samples <- list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/original", full.names=F) %>% gsub(".cov.gz","",.)
stopifnot(files.samples%in%sample_metadata$id_met)

for (i in files.samples) {
	
	tmp <- fread(grep(i,files,value=T), select=c(1,2,4)) %>%
		setnames(c("chr","pos","rate")) %>%
		.[rate!=50] %>%
		.[,rate:=round(rate/100)] %>%
		.[,chr:=paste0("chr",chr)]

	outfile <- file.path(outdir,sprintf("%s.tsv.gz",i))
	stopifnot(sprintf("%s.tsv.gz",i)%in%list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/met/cpg_level", full.names=F))
	fwrite(tmp, outfile, sep="\t", col.names=T, quote=F, na="NA")
}


##################
## scNMT-seq met ##
##################

indir <- "/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/met/cpg_level/old"
outdir <- "/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/met/cpg_level/new"
sample_metadata <- fread("/bi/group/reik/ricard/data/tet_chimera_nmtseq/results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")

files <- list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/met/cpg_level", full.names=T)

# Sanity check
files.samples <- list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/met/cpg_level", full.names=F, pattern="*.gz")
stopifnot(gsub(".tsv.gz","",files.samples)%in%sample_metadata$id_met)

for (i in files.samples) {
	tmp <- fread(file.path(indir,i), select=c(1,2,4)) %>% setnames(c("chr","pos","rate"))
	fwrite(tmp, file.path(outdir,i), sep="\t", col.names=T, quote=F, na="NA")
}


##################
## scMT-seq acc ##
##################

outdir <- "/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/met/cpg_level/new"
sample_metadata <- fread("/bi/group/reik/ricard/data/tet_chimera_nmtseq/results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")

files <- list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/acc/gpc_level", full.names=T, )

# Sanity check
files.samples <- list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/acc/gpc_level", full.names=F, pattern="*.gz") %>% gsub(".tsv.gz","",.)
samples.to.remove <- files.samples[!files.samples%in%sample_metadata$id_acc]

for (i in samples.to.remove) {
	file.remove(file.path(indir,sprintf("%s.tsv.gz",i)))
}


stopifnot(sample_metadata[!is.na(id_acc),id_acc] %in% gsub(".tsv.gz","",list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/acc/gpc_level", full.names=F, pattern="*.gz")))

##################
## scNMT-seq acc ##
##################

indir <- "/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/acc/gpc_level"
outdir <- "/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/acc/gpc_level/new"
sample_metadata <- fread("/bi/group/reik/ricard/data/tet_chimera_nmtseq/results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")

files <- list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/acc/cpg_level", full.names=T)

# Sanity check
files.samples <- list.files("/bi/group/reik/ricard/data/tet_chimera_nmtseq/processed/acc/gpc_level", full.names=F, pattern="*.gz")
stopifnot(gsub(".tsv.gz","",files.samples)%in%sample_metadata$id_acc)

for (i in files.samples) {
	tmp <- fread(file.path(indir,i), select=c(1,2,4)) %>% setnames(c("chr","pos","rate"))
	fwrite(tmp, file.path(outdir,i), sep="\t", col.names=T, quote=F, na="NA")
}
