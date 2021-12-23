here::i_am("rna/differential/differential.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$script <- here::here("rna/differential/differential.R")
io$outdir <- file.path(io$basedir,"results_new/rna/differential"); dir.create(io$outdir, showWarnings=F)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & !is.na(celltype.mapped) & ko_type!="crispr"] %>% 
  .[,ko:=ifelse(grepl("KO",sample),"TET TKO","WT")]

# Only consider cell types with sufficient observations in both KO and WT cells
celltypes.to.use <- sample_metadata[,.(N=.N),by=c("ko","celltype.mapped")] %>% .[N>=15] %>% .[,.N,by="celltype.mapped"] %>% .[N>1,celltype.mapped]
sample_metadata <- sample_metadata[celltype.mapped%in%celltypes.to.use]

print(table(sample_metadata$celltype.mapped,sample_metadata$ko))

###################################
## Run all pair-wise comparisons ##
###################################

for (i in celltypes.to.use) {
  outfile <- sprintf("%s/%s_WT_vs_KO.txt.gz", io$outdir,i)
  
  if (!file.exists(outfile)) {
    
    # Define LSF command
    if (grepl("BI",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
      lsf <- sprintf("sbatch -n 1 --mem 8G --wrap")
    }
    cmd <- sprintf("%s 'Rscript %s --celltype %s --groupA WT --groupB KO --group_label ko --outfile %s'", lsf, io$script, i, outfile)
    
    # Run
    print(cmd)
    system(cmd)
  }
}

