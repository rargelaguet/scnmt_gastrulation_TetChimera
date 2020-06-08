
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")
}  
io$mapping.dir <- paste0(io$basedir,"/results/iterative_mapping")
io$outdir <- paste0(io$basedir,"/results/iterative_mapping/pdf")

opts$batches <- c(
  # E7.5
  # "E75_TET_TKO_L002",
  # "E75_WT_Host_L001",
  
  # E8.5
  "E85_Rep1_TET_TKO_L004",
  "E85_Rep2_TET_TKO_L006",
  "E85_Rep1_WT_Host_L003",
  "E85_Rep2_WT_Host_L005"
  
  # E12.5
  # "E125_DNMT3A_HET_A_L001",
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002",
  # "E125_DNMT3A_KO_E_L004"
)

###################
## Load metadata ##
###################

sample_metadata <- sample_metadata %>%
  .[batch%in%opts$batches] %>%
table(sample_metadata$batch)

################
## Load data  ##
################

mapping_dt <- opts$batches %>% map(function(i) {
  rbind(
    fread(sprintf("%s/%s_standard_mnn.txt.gz",io$mapping.dir,i)) %>% .[,batch:=i] %>% .[,method:="Standard MNN"],
    fread(sprintf("%s/%s_iterative_mnn.txt.gz",io$mapping.dir,i)) %>% .[,batch:=i] %>% .[,method:="Tree-guided MNN"]
  )
}) %>% rbindlist

##########
## Plot ##
##########

barplot.pub <- function(df, x, colors=NULL, xlim.max=NULL) {
  p <- ggplot(df, aes_string(x=x, y="N")) +
    scale_x_discrete(drop=FALSE) + 
    labs(y="Number of cells") +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.3)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.3), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  if (is.null(colors)) {
    p <- p + geom_bar(stat="identity", color="black")
  } else {
    p <- p + geom_bar(aes_string(fill=x), stat="identity", color="black") + 
      scale_fill_manual(values=colors, drop=F)
  }
  
  if (!is.null(xlim.max)) {
    p <- p + coord_flip(ylim=c(0,xlim.max))
  } else {
    p <- p + coord_flip()
  }
  
  return(p)
}

for (i in opts$batches) {
  
  to.plot <- mapping_dt[batch==i] %>% 
    merge(sample_metadata, by=c("cell","batch")) %>% 
    .[celltype_mapped=="Forebrain Midbrain Hindbrain",celltype_mapped:="Forebrain_Midbrain_Hindbrain"] %>%
    .[,.N,by=c("celltype_mapped","batch","method")]
  
  to.plot[,celltype_mapped:=stringr::str_replace_all(celltype_mapped," ", "_")]
  
  p <- barplot.pub(to.plot, x="celltype_mapped", colors=opts$celltype.colors) +
    facet_wrap(~method, nrow=1, scales="free_x") +
    theme(
      strip.background = element_blank()
    )  
  
  print(p)
}
