here::here("metacc/quantify_feature_level/quantify_feature_level.R")

suppressMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='Cell metadata')
p$add_argument('--indir',  type="character",              help='Input directory')
p$add_argument('--outdir',  type="character",              help='Output directory')
p$add_argument('--featuresdir',  type="character",              help='Features directory')
p$add_argument('--context',  type="character",              help='cg/CG or gc/GC')
p$add_argument('--annos',    type="character",  nargs="+",  help='Genomic annotation')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

###################
## Load settings ##
###################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
args <- list()
args$indir <- file.path(io$basedir,"processed/met/cpg_level")
args$outdir <- file.path(io$basedir,"processed/met/feature_level")
args$features.dir  <- file.path(io$basedir,"features/filt")
args$metadata <- file.path(io$basedir,"results/met/qc/sample_metadata_after_metQC.txt.gz")
args$context <- "CG"
args$annos <- c("prom_2000_2000")
## END TEST ##

# Sanity checks
stopifnot(args$context %in% c("CG","GC"))

#############
## Options ##
#############

# Define cells
if (args$context=="CG") {
  samples_keep <- fread(args$metadata) %>% .[!is.na(id_met),id_met]
} else if (args$context=="GC") {
  samples_keep <- fread(args$metadata) %>% .[!is.na(id_acc),id_acc]
}

# ###########
# ## Print ##
# ###########

cat("\nProcessing methylation samples with the following options:\n")
cat(sprintf("- Input folder for annotation: %s\n",args$features.dir))
cat(sprintf("- Input folder for bismark files: %s\n",args$indir))
cat(sprintf("- Output folder: %s\n",args$outdir))
cat(sprintf("- Annotations: %s\n", paste(args$annos, collapse=" ")))
cat(sprintf("- Annotating CG or GC?: %s\n",args$context))

######################
## Load annotations ##
######################

anno_list <- list()
for (i in 1:length(args$annos)) {
  anno_list[[i]] <- fread(file.path(args$features.dir,"%s.bed.gz",args$annos[i]), header=F, select=c(1,2,3,4,5)) %>%
    setnames(c("chr","start","end","strand","id")) %>%
    .[,chr:=as.factor(chr)] %>%
    setkey(chr,start,end)
}
names(anno_list) <- args$annos


#########################################
## Preprocess and annotate the samples ##
#########################################

# Create ouput temporary folder
dir.create(sprintf("%s/tmp",args$outdir), recursive=T)

for (i in 1:length(samples_keep)) {
  files_processed <- list.files(sprintf("%s/tmp",args$outdir))
  if (all(sprintf("%s_%s.gz",samples_keep[i],args$annos) %in% files_processed)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",samples_keep[i])) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",samples_keep[i]))  
    
    # Read and parse raw methylation data
    dat_sample <- fread(sprintf("%s/%s.tsv.gz",args$indir,samples_keep[i]), sep="\t", verbose=F, showProgress=F) %>%
      .[,c("chr","pos","rate")] %>%
      .[,c("start","end") := list(pos,pos)] %>% # Add 'start' and 'end' columns to do the overlap
      .[,c("chr","pos"):=list(as.factor(chr),NULL)] %>% 
      setkey(chr,start,end)
    
    # Sanity check
    stopifnot(all(dat_sample$rate %in% c(0,1)))
    
    # Overlap data with annotations
    for (anno in args$annos) {
      fname.out <- sprintf("%s/tmp/%s_%s.gz",args$outdir,samples_keep[i],anno)
      if (file.exists(fname.out)) {
        cat(sprintf("Annotation for %s with %s already found, loading...\n",samples_keep[i],anno))
      } else {
        cat(sprintf("Annotating %s with %s annotations...\n",samples_keep[i],anno))
        
        # Overlap and calculate methylation status for each region in the annotation by summarising over all sites
        ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch=0) %>% 
          .[,"i.end":=NULL] %>% setnames("i.start","pos") %>%
          .[,c("sample","anno") := list(samples_keep[i],anno)] %>%
          # Compute number of methylated CpGs and the corresponding methylation rates
          .[,.(rate=round(mean(rate)*100), Nmet=sum(rate==1), N=.N), keyby=c("sample","id","anno")] %>%
          # Reorder columns
          .[,c("sample","id","anno","Nmet","N","rate")]
        
        # Store and save results
        fwrite(ov, fname.out, quote=FALSE, sep="\t", col.names=FALSE)
      }
    }
  }
}

#####################################
## Concatenate everything and save ##
#####################################

for (i in args$annos) {
  outfile <- sprintf("%s/%s.tsv.gz", args$outdir, i)
  if(file.exists(outfile)) {
    cat(sprintf("File %s already exists, ignoring...\n", outfile))
  } else {
    system(sprintf("cat %s/tmp/*_%s.gz | zgrep -E -v 'sample|id_met|id_acc' | pigz > %s", args$outdir, i, outfile))
  }
}
