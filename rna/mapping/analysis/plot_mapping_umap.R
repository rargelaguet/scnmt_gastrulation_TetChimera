here::i_am("rna/mapping/analysis/plot_mapping_umap.R")

source(here::here("settings.R"))
source(here::here("rna/mapping/analysis/plot_utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--query_metadata',        type="character",                               help='Cell metadata (after mapping)')
# p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--atlas_metadata',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$query_metadata <- file.path(io$basedir,"results_new/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$atlas_metadata <- file.path(io$atlas.basedir,"sample_metadata.txt.gz")
args$outdir <- file.path(io$basedir,"results_new/rna/mapping/pdf")
## END TEST ##

dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"per_sample"), showWarnings = F)
dir.create(file.path(args$outdir,"per_class"), showWarnings = F)

#####################
## Define settings ##
#####################

# Options

# Dot size
opts$size.mapped <- 0.30
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.65
opts$alpha.nomapped <- 0.35

opts$remove_ExE_cells <- TRUE
opts$subset_atlas_cells <- TRUE

#########################
## Load query metadata ##
#########################

sample_metadata <- fread(args$query_metadata) %>%
  .[,class2:=ifelse(grepl("WT",class),"WT","Tet-TKO")] %>%
  # .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(closest.cell)]
  .[pass_rnaQC==TRUE & !is.na(closest.cell)]

stopifnot("closest.cell"%in%colnames(sample_metadata))

if (opts$remove_ExE_cells) {
  sample_metadata <- sample_metadata %>% .[!celltype.mapped%in%c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

################
## Load atlas ##
################

# Load atlas cell metadata
meta_atlas <- fread(args$atlas_metadata) %>%
  # .[celltype%in%opts$celltypes] %>%
  .[stripped==F & doublet==F]

if (opts$remove_ExE_cells) {
  meta_atlas <- meta_atlas %>% .[!celltype%in%c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

# Subset cells to speed up plotting
if (opts$subset_atlas_cells) {
  meta_atlas <- rbind(
    meta_atlas[cell%in%unique(sample_metadata$closest.cell)],
    meta_atlas %>% .[!cell%in%unique(sample_metadata$closest.cell)] %>% .[sample.int(n=nrow(.), size=nrow(.)/4)]
  )
}

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas %>%
  .[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

####################
## Plot all cells ##
####################

# to.plot <- umap.dt %>% copy %>%
#   .[,index:=match(cell, sample_metadata[,closest.cell] )] %>% 
#   .[,mapped:=as.factor(!is.na(index))] %>% 
#   .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas cells","Query cells"))] %>%
#   setorder(mapped) 
# 
# p <- plot.dimred(to.plot, query.label = "Query cells", atlas.label = "Atlas cells")
# 
# pdf(sprintf("%s/umap_mapped_allcells.pdf",args$outdir), width=8, height=6.5)
# print(p)
# dev.off()

###############################
## Plot one sample at a time ##
###############################

samples.to.plot <- unique(sample_metadata$sample)

for (i in samples.to.plot) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[sample==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/per_sample/umap_mapped_%s.pdf",args$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

###############################
## Plot one class at a time ##
###############################

classes.to.plot <- unique(sample_metadata$class)

for (i in classes.to.plot) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[class==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas") + theme(legend.position = "none")
  
  pdf(sprintf("%s/per_class/umap_mapped_%s.pdf",args$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

###########################################
## Plot multiple clases at the same time ##
###########################################

samples.to.plot <- c(
  "E7.5_WT",
  "E7.5_TET_TKO",
  # "E7.5_TET_TKO_crispr",
  "E8.5_WT",
  "E8.5_WT_CD41+",
  "E8.5_WT_KDR+",
  "E8.5_WT_KDR+_CD41+",
  "E8.5_TET_TKO",
  "E8.5_TET_TKO_KDR+",
  "E8.5_TET_TKO_CD41+",
  "E8.5_TET_TKO_KDR+_CD41+"
)
# samples.to.plot <- unique(sample_metadata$sample)

to.plot <- samples.to.plot %>% map(function(i) {
  umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[sample==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas","Query"))] %>%
    .[,sample:=i] %>%
    setorder(mapped)
}) %>% rbindlist %>% .[,sample:=factor(sample,levels=samples.to.plot)]

p <- plot.dimred(to.plot, query.label = "Query", atlas.label = "Atlas") + 
  facet_wrap(~sample, nrow=2) +
  theme(legend.position = "none") 

pdf(file.path(args$outdir,"per_class/umap_mapped_all_classes.pdf"), width=12, height=6)
print(p)
dev.off()

#############################
## Plot WT and KO together ##
#############################

# Subsample query cells to have the same N per class
# sample_metadata_subset <- sample_metadata %>% .[,.SD[sample.int(n=.N, size=4500)], by=c("stage","class2")]

# i <- "E7.5"
to.plot <- umap.dt %>% copy %>%
  .[,index.wt:=match(cell, sample_metadata[class2=="WT",closest.cell] )] %>%
  .[,index.ko:=match(cell, sample_metadata[class2=="Tet-TKO",closest.cell] )] %>%
  .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
  .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
  .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
  .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas","WT","Tet-TKO"))] %>% setorder(mapped)

p <- plot.dimred.wtko(to.plot, wt.label = "WT", ko.label = "Tet-TKO", nomapped.label = "Atlas") +
  theme(legend.position = "top", axis.line = element_blank())

pdf(sprintf("%s/per_class/umap_mapped_WT_and_KO.pdf",args$outdir), width=5.5, height=6.5)
print(p)
dev.off()

# Completion token
file.create(file.path(args$outdir,"completed.txt"))
