suppressPackageStartupMessages(library(argparse))

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

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$samples <- opts$samples
# args$celltype_label <- "celltype.mapped"
# args$outdir <- file.path(io$basedir,"results/rna/celltype_proportions")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)

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

################################################
## Calculate cell type proportions per sample ##
################################################

to.plot <- sample_metadata %>%
  .[,N:=.N,by="sample"] %>%
  # .[,celltype:=stringr::str_replace_all(celltype,opts$aggregate.celltypes)] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","celltype")] %>%
  setorder(sample)  %>% .[,sample:=factor(sample,levels=args$samples)]

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

#########################
## Horizontal barplots ##
#########################

# Rename "_" to " " in cell types
# to.plot[,celltype:=stringr::str_replace_all(celltype,"_"," ")]
# names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=names(opts$celltype.colors))]


for (i in unique(to.plot$sample)) {
  p <- ggplot(to.plot[sample==i], aes(x=celltype, y=N)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    scale_x_discrete(drop=F) +
    # facet_wrap(~sample, nrow=1, scales="fixed") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.5)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    )
  
  pdf(file.path(args$outdir,sprintf("celltype_proportions_%s_horizontal_barplots.pdf",i)), width=7, height=5)
  print(p)
  dev.off()
}
