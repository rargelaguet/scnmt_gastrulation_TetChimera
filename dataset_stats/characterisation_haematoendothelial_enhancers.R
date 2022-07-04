# here::i_am("")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"features/other"); dir.create(io$outdir, showWarnings=F)

# Options
opts$celltypes.atlas.to.plot <- c(
  "Haematoendothelial_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3"
)

opts$rename.celltypes <- c(
  "Erythroid1" = "early_Erythroid",
  "Erythroid2" = "early_Erythroid",
  "Erythroid3" = "late_Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
)

opts$celltypes.scnmt.to.plot <- c(
  "Haematoendothelial_progenitors",
  "Blood_progenitors",
  "early_Erythroid",
  "late_Erythroid"
)

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$atlas.marker_genes) %>% 
  .[celltype%in%opts$celltypes.atlas.to.plot] %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename.celltypes)] %>%
  .[,.(logFC=mean(logFC), score=mean(score)), by=c("celltype","gene")] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes.scnmt.to.plot)]

# Merge cell types
table(marker_genes.dt$celltype)

#######################
## Load marker peaks ##
#######################

marker_peaks.dt <- fread(io$multiome.marker_peaks) %>% 
  .[celltype%in%opts$celltypes.atlas.to.plot] %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename.celltypes)] %>%
  .[,.(logFC=mean(logFC), score=mean(score)), by=c("celltype","feature")] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes.scnmt.to.plot)]

table(marker_peaks.dt$celltype)

########################
## Load peak metadata ##
########################

peak_metadata.dt <- fread(io$multiome.peak_metadata) %>%
  .[,peak:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[peak%in%unique(marker_peaks.dt$feature)] %>%
  .[,peakType:=factor(peakType,levels=c("Promoter","Intronic","Exonic","Distal"))]

table(peak_metadata.dt$peakType)

##########################
## Load peak2gene links ##
##########################

peak2genes.dt <- fread(io$multiome.peak2genes.all) %>% 
  .[peak%in%unique(marker_peaks.dt$feature) & gene%in%unique(marker_genes.dt$gene)]

length(unique(peak2genes.dt$gene))
length(unique(peak2genes.dt$peak))

#####################
## Load DE results ##
#####################

io$diff_results.folder <- "/Users/rargelaguet/data/10x_gastrulation_TetChimera/results_all/differential"
# opts$celltypes.scnmt.to.plot
diff_results.dt <- fread(file.path(io$diff_results.folder,"Blood_progenitors_WT_vs_TET_TKO.txt.gz")) %>%
  .[gene%in%marker_genes.dt$gene]

# Erythroid_WT_vs_TET_TKO.txt.gz
# Haematoendothelial_progenitors_WT_vs_TET_TKO.txt.gz
# Blood_progenitors_WT_vs_TET_TKO.txt.gz

###############################################
## Plot number of marker peaks per cell type ##
###############################################

to.plot <- marker_peaks.dt %>% .[,.N,by="celltype"] 

p <- ggbarplot(to.plot, x="celltype", y="N", fill="celltype") +
  labs(x="", y="Number of marker peaks") +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(colour="black",size=rel(0.8)),  
    # axis.text.x = element_text(colour="black",size=rel(0.65), angle=45, hjust=1, vjust=1),  
    axis.text.x = element_text(colour="black",size=rel(0.65))
  )

pdf(file.path(io$outdir,"barplot_number_marker_peaks_per_celltype.pdf"), width = 6, height = 4)
print(p)
dev.off()

#####################################
## Plot fraction of peaks per type ##
#####################################

to.plot <- peak_metadata.dt %>%
  .[,.N,by="peakType"] %>% .[,percentage:=100*N/sum(N)]

p <- ggpie(to.plot, x="N", label="peakType", fill="peakType") +
  theme(
    legend.position = "none"
  )

pdf(file.path(io$outdir,"pieplot_atac_peak_type.pdf"), width = 5, height = 4)
print(p)
dev.off()

############################################################
## Plot number of associations between peaks and DE genes ##
############################################################

tmp <- peak2genes.dt[,c("gene","peak","dist")] %>% 
  merge(diff_results.dt[,c("gene","logFC","padj_fdr")], by="gene", all=TRUE)

mean(tmp[,any(padj_fdr<=0.01,na.rm=T),by="peak"][[2]])
mean(tmp[,any(padj_fdr<=0.01,na.rm=T),by="gene"][[2]])
