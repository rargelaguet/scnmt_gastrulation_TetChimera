##############
## Settings ##
##############

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
source("/Users/ricard/10x_gastrulation_TetChimera/differential/analysis/utils.R")
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/differential/pdf")

# Define groups
opts$groupA <- c("E8.5_Host")
opts$groupB <- c("E8.5_TET_TKO")

#############################################
## Load results from differential analysis ##
#############################################

dt <- opts$celltypes %>% map(function(i) {
  file <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$diff.dir,i,opts$groupA,opts$groupB)
  if (file.exists(file)) fread(file) %>% .[,c("celltype","groupA","groupB"):=list(i,opts$groupA,opts$groupB)]
}) %>% rbindlist

# (TO-DO) Filter by minimum number of cells per group

for (i in unique(dt$celltype)) {
  to.plot <- dt[celltype==i] %>% .[!is.na(sig)]
  p <- gg_volcano_plot(to.plot, top_genes=10)
  
  pdf(sprintf("%s/%s_%s_%s_volcano.pdf",io$outdir,i,opts$groupA,opts$groupB), width=9, height=5, useDingbats = F)
  print(p)
  dev.off()
}
