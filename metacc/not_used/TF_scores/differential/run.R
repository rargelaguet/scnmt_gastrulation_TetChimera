#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))

io$script <- here::here("metacc/differential/feature_level/diff_metacc_feature_level.R")
# io$tmpdir <- file.path(io$basedir,"results_new/atac/archR/differential/tmp"); dir.create(io$tmpdir, showWarnings=F)

opts$group_label <- "ko"

opts$min_cells_total <- 15
opts$min_cells_per_feature <- 10

# Genomic contexts
# opts$annos <- c("multiome_peaks")
opts$annos <- NULL
if (is.null(opts$annos)) {
  opts$annos <- list.files(io$met_data_parsed, pattern=".tsv.gz") %>% gsub(".tsv.gz","",.)# %>% head(n=1)
}

opts$celltypes <- c(
  "Primitive_Streak",
  "Nascent_mesoderm",
  "Mesenchyme",
  "ExE_mesoderm",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors",
  # "Surface_ectoderm",
  "early_Erythroid"
)

opts$group_label <- "ko"

#########
## Run ##
#########


io$outdir <- c(
  "CG" = file.path(io$basedir,"results_new/met/differential"),
  "GC" = file.path(io$basedir,"results_new/acc/differential")
)


groupA <- "WT"; groupB <- "KO"

for (context in c("CG","GC")) {
  
  if (context=="CG") {
    sample_metadata <- fread(io$metadata) %>%  .[pass_metQC==TRUE & !is.na(celltype.mapped)] %>% .[,ko:=ifelse(grepl("KO",class),"KO","WT")]
  } else {
    sample_metadata <- fread(io$metadata) %>%  .[pass_accQC==TRUE & !is.na(celltype.mapped)] %>% .[,ko:=ifelse(grepl("KO",class),"KO","WT")]
  }
  
  celltypes <- sample_metadata[,.N,by=c("celltype.mapped","ko")] %>% .[,.(N=sum(N>=opts$min_cells_total)),by="celltype.mapped"] %>% .[N==2,celltype.mapped]
  for (i in celltypes) {
    for (anno in opts$annos) {
      outfile <- file.path(io$outdir[[context]],sprintf("%s_%s_%s_vs_%s.txt.gz", anno, i, groupA, groupB))
      if (file.exists(outfile)) {
        print(sprintf("%s already exists, skipping...",outfile))
      } else {
  
        # Define LSF command
        if (grepl("BI",Sys.info()['nodename'])) {
          lsf <- ""
        } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
          lsf <- sprintf("sbatch -n 1 --mem 9G --wrap")
        }
        cmd <- sprintf("%s 'Rscript %s --groupA %s --groupB %s --celltype %s --anno %s --group_label %s --context %s --min_cells %d --outfile %s'", 
          lsf, io$script, groupA, groupB, i, anno, opts$group_label, context, opts$min_cells_per_feature, outfile)
  
        # Run
        print(cmd)
        system(cmd)
        
      }
    }
  cat("\n\n")
  }
}
