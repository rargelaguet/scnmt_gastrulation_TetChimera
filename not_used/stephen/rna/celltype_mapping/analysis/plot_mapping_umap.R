#####################
## Define settings ##
#####################

# load default setings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna/mapping/analysis/plot_utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/rna/mapping/analysis/plot_utils.R")
} else {
  source(here::here("settings.R"))
  source(here::here("rna/celltype_mapping/analysis/plot_utils.R"))
}

# I/O
io$outdir <- paste0(io$basedir,"/results/rna/mapping/pdf")
dir.create(io$outdir, recursive = TRUE, showWarnings = F)
dir.create(paste0(io$outdir,"/per_sample"), showWarnings = F)
dir.create(paste0(io$outdir,"/per_stage"), showWarnings = F)

# Options

# Dot size
opts$size.mapped <- 0.2
opts$size.nomapped <- 0.01

# Transparency
opts$alpha.mapped <- 1
opts$alpha.nomapped <- 0.1

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE] #%>%
 # .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
meta_atlas <- fread(io$atlas.metadata) %>%
  .[celltype%in%opts$celltypes] %>%
  .[stripped==F & doublet==F]

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

#########################################################
## Plot dimensionality reduction: one sample at a time ##
#########################################################
plates <- sample_metadata[, unique(plate)]
plot_list <- vector("list", length = length(plates))
names(plot_list) <- plates

for (i in plates) {
  
  meta <- sample_metadata[plate == i]
  
  meta_types <- meta[, .(ko_type, stage, tdTOM, `KDR-Cy7`, `CD41-BV421`)] %>%
    unique()
  
  sample_type <- paste(
    meta_types$stage, 
    meta_types$ko_type, 
    ifelse(meta_types$tdTOM, "tdTom+", "tdTom-"),
    ifelse(meta_types$`KDR-Cy7`, "KDR+", "KDR-"),
    ifelse(meta_types$`CD41-BV421`, "CD41+", "CD41-"),
    sep = " "
    ) 
  
  # not working for some reason???
  
  # to.plot <- umap.dt %>% copy %>%
  #   .[,index:=match(cell, sample_metadata[plate==i,closest.cell] )] %>% 
  #   .[,mapped:=as.factor(!is.na(index))] %>% 
  #   .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
  #   setorder(mapped) 
  
  to.plot <- copy(umap.dt) %>%
    merge(meta[, .(cell = closest.cell, mapped = sample_type)], by = "cell", all.x = TRUE) %>%
    .[is.na(mapped), mapped := "Atlas"] 
  
  p <- plot.dimred(to.plot, query.label = sample_type, atlas.label = "Atlas")
  
  pdf(sprintf("%s/per_sample/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
  
  plot_list[[i]] <- p
}


big_plot <- cowplot::plot_grid(plotlist = plot_list)
big_plot
########################################################
## Plot dimensionality reduction: one stage at a time ##
########################################################

for (i in opts$stages) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[stage==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/per_stage/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}
