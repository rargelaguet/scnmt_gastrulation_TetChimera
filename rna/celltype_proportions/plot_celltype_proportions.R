here::i_am("rna/celltype_proportions/plot_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--celltype_label', type="character", help='Cell type label')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################


## START TEST ##
# args$metadata <- file.path(io$basedir,"results_new/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$samples <- opts$samples
# args$celltype_label <- "celltype.mapped"
# args$outdir <- file.path(io$basedir,"results_new/rna/celltype_proportions")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"per_sample"), showWarnings = F)
dir.create(file.path(args$outdir,"per_plate"), showWarnings = F)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%args$samples & !is.na(eval(as.name(args$celltype_label)))] %>%
  setnames(args$celltype_label,"celltype")

# if (!"stage"%in%colnames(sample_metadata)) {
#   sample_metadata[,stage:=stringr::str_replace_all(sample,opts$sample2stage)]
# }
# sample_metadata[,stage:=stringr::str_replace_all(sample,opts$sample2stage)]

table(sample_metadata$sample)
table(sample_metadata$plate)

######################
## Stacked barplots ##
######################

# p <- ggplot(to.plot, aes(x=sample, y=celltype_proportion)) +
#   geom_bar(aes(fill=celltype), stat="identity", color="black") +
#   facet_wrap(~class, scales = "free_x", nrow=1) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.title = element_blank(),
#     axis.text.y = element_blank(),
#     # axis.text.x = element_text(color="black", size=rel(0.75)),
#     axis.text.x = element_blank(),
#     axis.ticks = element_blank(),
#     axis.line = element_blank()
#   )
# 
# pdf(file.path(args$outdir,"celltype_proportions_stacked_barplots.pdf"), width=7, height=5)
# print(p)
# dev.off()

########################
## Barplot per sample ##
########################

# Calculate cell type proportions per sample
celltype_proportions.dt <- sample_metadata %>%
  .[,N:=.N,by="sample"] %>%
  # .[,celltype:=stringr::str_replace_all(celltype,opts$aggregate.celltypes)] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","celltype")] %>%
  setorder(sample)  %>% .[,sample:=factor(sample,levels=opts$samples)]

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(celltype_proportions.dt$celltype)]
celltype_proportions.dt[,celltype:=factor(celltype, levels=rev(names(opts$celltype.colors)))]

samples.to.plot <- unique(celltype_proportions.dt$sample)

for (i in samples.to.plot) {
  
  to.plot <- celltype_proportions.dt[sample==i]
  
  p <- ggplot(to.plot, aes(x=celltype, y=N)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    # facet_wrap(~plate, nrow=1, scales="fixed") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.9)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    )
  
  pdf(sprintf("%s/per_sample/celltype_proportions_%s.pdf",args$outdir,i), width=7, height=5)
  print(p)
  dev.off()
}

########################
## Barplots per plate ##
########################

celltype_proportions.dt <- sample_metadata %>%
  .[,N:=.N,by=c("sample","plate")] %>%
  # .[,celltype:=stringr::str_replace_all(celltype,opts$aggregate.celltypes)] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","celltype","plate")] %>%
  setorder(sample)  %>% .[,sample:=factor(sample,levels=opts$samples)]

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(celltype_proportions.dt$celltype)]
celltype_proportions.dt[,celltype:=factor(celltype, levels=rev(names(opts$celltype.colors)))]

samples.to.plot <- unique(celltype_proportions.dt$sample)

for (i in samples.to.plot) {
  
  to.plot <- celltype_proportions.dt[sample==i] %>% 
    .[,N:=sum(N),by="celltype"] %>% .[N>=5]
  
  p <- ggplot(to.plot[sample==i], aes(x=celltype, y=N)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~plate, scales="fixed") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.75)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    )
  
  pdf(sprintf("%s/per_plate/celltype_proportions_%s_per_plate.pdf",args$outdir,i), width=7, height=7)
  print(p)
  dev.off()
}
