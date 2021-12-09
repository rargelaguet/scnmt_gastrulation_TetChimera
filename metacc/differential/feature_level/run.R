#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))

io$script <- here::here("metacc/differential/feature_level/diffmet_feature_level.R")
# io$tmpdir <- file.path(io$basedir,"results_new/atac/archR/differential/tmp"); dir.create(io$tmpdir, showWarnings=F)

opts$group_label <- "celltype_class"

opts$min_cells <- 10


# Genomic contexts
opts$annos <- c("prom_2000_2000")
# opts$annos <- NULL
# if (is.null(opts$annos)) {
#   opts$annos <- list.files(io$met_data_parsed, pattern=".tsv.gz") %>% gsub(".tsv.gz","",.)# %>% head(n=3)
# }

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[,celltype_class:=sprintf("%s_%s",celltype.mapped_mnn,class)] %>% # temporary
  .[pass_metQC==TRUE]

stopifnot(opts$group_label%in%colnames(sample_metadata))

sample_metadata <- sample_metadata %>%
  setnames(opts$group_label,"group")

# subset groups with sufficient number of cells
sample_metadata <- sample_metadata %>%
  .[,N:=.N,by="group"] %>% .[N>=opts$min_cells] %>% .[,N:=NULL]

table(sample_metadata$group)

#########
## Run ##
#########

opts$celltypes <- unique(sample_metadata$celltype.mapped_mnn) %>% head(n=2)

io$outdir <- file.path(io$basedir,sprintf("results/met/differential/%s",opts$group_label))

# i <- "Blood_progenitors_2"; anno <- "prom_2000_2000"
for (i in opts$celltypes) {
  tmp <- sample_metadata[celltype.mapped_mnn==i,.N,by="group"]
  if (all(tmp$N>=opts$min_cells)) {
    groupA <- tmp$group[[1]]; groupB <- tmp$group[[2]]
    
    for (anno in opts$annos) {
      
      outfile <- file.path(io$outdir,sprintf("%s_%s_vs_%s.txt.gz", anno, groupA, groupB))
      if (file.exists(outfile)) {
        print(sprintf("%s already exists, skipping...",outfile))
      } else {
  
        # Define LSF command
        if (grepl("BI",Sys.info()['nodename'])) {
          lsf <- ""
        } else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
          lsf <- sprintf("sbatch -n 1 --mem 7G --wrap")
        }
        cmd <- sprintf("%s 'Rscript %s --groupA %s --groupB %s --anno %s --group_label %s --context CG --min_cells %d --outfile %s'", 
          lsf, io$script, groupA, groupB, anno, opts$group_label, opts$min_cells, outfile)
  
        # Run
        print(cmd)
        system(cmd)
        
      }
    }
  }
  
}
