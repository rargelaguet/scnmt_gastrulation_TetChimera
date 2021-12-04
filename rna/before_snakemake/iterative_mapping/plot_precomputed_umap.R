
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_spatial/settings.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation_spatial/settings.R")
}  
io$mapping.dir <- paste0(io$basedir,"/iterative_mapping")
io$outdir <- paste0(io$basedir,"/iterative_mapping/pdf")

opts$embryos <- c("embryo1","embryo2","embryo3")

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[embryo%in%opts$embryos] 

################
## Load data  ##
################

mapping_dt <- opts$embryos %>% map(function(i) {
  rbind(
    fread(sprintf("%s/%s_standard_mnn.txt.gz",io$mapping.dir,i)) %>% .[,embryo:=i] %>% .[,method:="Standard MNN"],
    fread(sprintf("%s/%s_iterative_mnn.txt.gz",io$mapping.dir,i)) %>% .[,embryo:=i] %>% .[,method:="Tree-guided MNN"]
  )
}) %>% rbindlist

##########
## Plot ##
##########

for (i in opts$embryos) {
  to.plot <- mapping_dt[embryo==i] %>% 
    merge(sample_metadata[,c("cell","embryo","z","x_global_affine","y_global_affine")], by=c("cell","embryo")) %>% 
    setnames(c("x_global_affine","y_global_affine"),c("V1","V2")) %>%
    .[celltype_mapped=="Forebrain Midbrain Hindbrain",celltype_mapped:="Forebrain/Midbrain/Hindbrain"]# %>%
    # setnames(c("Tree-guided MNN","MNN")) %>%
    # melt(id.vars=c("cell","embryo","z","V1","V2"), measure.vars=c("Tree-guided MNN","MNN"), variable.name="class", value.name="celltype")
  
  p <- ggplot(to.plot, aes(x=V1, y=-V2)) +
    geom_point(aes(fill=celltype_mapped), size=1.25, shape=21, color="black", stroke=0.05) +
    facet_grid(z~method) +
    scale_fill_manual(values=opts$celltype.colors) +
    guides(fill = guide_legend(override.aes = list(size=4))) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
  
  pdf(sprintf("%s/umap_%s.pdf",io$outdir,i), width=14, height=9, useDingbats = F)
  print(p)
  dev.off()
}
