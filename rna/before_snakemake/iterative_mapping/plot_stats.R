
#####################
## Define settings ##
#####################

source("/Users/ricard/scnmt_gastrulation_TetChimera/settings.R")
source("/Users/ricard/scnmt_gastrulation_TetChimera/rna/iterative_mapping/plot_utils.R")
io$mapping.dir <- paste0(io$basedir,"/rna/results/iterative_mapping")
io$outdir <- paste0(io$basedir,"/rna/results/iterative_mapping/pdf")

opts$plates <- c(
  # "tet_chimera_march20_plate1"
  # "tet_chimera_march20_plate2"
  # "tet_chimera_march20_plate3"
  # "tet_chimera_march20_plate4"
  # "tet_chimera_march20_plate5",
  # "tet_chimera_march20_plate6"
    "tet_chimera_oct19_plate1",
    "tet_chimera_oct19_plate2",
    "tet_chimera_oct19_plate4",
    "tet_chimera_oct19_plate5"
)


################### 
## Load metadata ##
###################

# sample_metadata <- sample_metadata %>%
sample_metadata <- fread("/Users/ricard/data/scnmt_gastrulation_TetChimera/TO_MERGE_scnmt_gastrulation_TetKO/sample_metadata.txt.gz") %>%
  .[plate%in%opts$plates] %>%
  setnames("id_rna","cell")
table(sample_metadata$plate)

################
## Load data  ##
################

mapping_dt <- opts$plates %>% map(function(i) {
  rbind(
    # fread(sprintf("%s/%s_standard_mnn.txt.gz",io$mapping.dir,i)) %>% .[,plate:=i] %>% .[,method:="Standard MNN"],
    fread(sprintf("%s/%s_iterative_mnn.txt.gz",io$mapping.dir,i)) %>% .[,plate:=i] %>% .[,method:="Tree-guided MNN"]
  )
}) %>% rbindlist

##########
## Plot ##
##########

p_list <- list()
for (i in opts$plates) {
  
  to.plot <- mapping_dt[plate==i] %>% 
    merge(sample_metadata, by=c("cell","plate")) %>% 
    .[celltype_mapped=="Forebrain Midbrain Hindbrain",celltype_mapped:="Forebrain_Midbrain_Hindbrain"] %>%
    .[,.N,by=c("celltype_mapped","plate","method","tdTOM")]
  
  to.plot[,celltype_mapped:=stringr::str_replace_all(celltype_mapped," ", "_")]
  
  p_list[[i]] <- barplot.pub(to.plot, x="celltype_mapped", colors=opts$celltype.colors) +
    ggtitle(i) +
    # facet_wrap(~tdTOM, nrow=1, scales="free_x") +
    facet_wrap(~tdTOM, nrow=2, scales="free_x") +
    theme(
      strip.background = element_blank()
    )  
  
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=8, height=6.5, useDingbats = F)
  # print(p_list[[i]])
  # dev.off()
}

cowplot::plot_grid(plotlist=p_list)

