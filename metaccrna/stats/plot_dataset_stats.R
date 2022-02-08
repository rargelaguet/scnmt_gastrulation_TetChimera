source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results_new/metaccrna/stats")
dir.create(io$outdir, showWarnings=F)

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[,ko:=ifelse(grepl("KO",class),"Tet-TKO","WT")] %>%
  .[,ko:=factor(ko,levels=c("WT","Tet-TKO"))]

###################################################
## Number of cells with multi-omics measurements ##
###################################################

to.plot <- sample_metadata %>%
  .[,.(met=sum(!is.na(id_met)), acc=sum(!is.na(id_acc)), rna=sum(!is.na(id_rna))),by="ko"] %>%
  melt(id.vars="ko")

facet.labels <- c(rna = "RNA expression", met = "DNA methylation", acc = "Chr. accessibility")
  
p <- ggplot(to.plot, aes_string(x="ko", y="value", fill="ko")) +
  geom_bar(stat="identity", position="dodge", color="black") +
  facet_wrap(~variable, labeller = as_labeller(facet.labels)) +
  labs(x="", y="Number of cells") +
  # scale_fill_manual(values=c(acc="#00BFC4", rna="#00CD00", met="#F8766D")) +
  scale_fill_manual(values=opts$class.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.y = element_text(size=rel(1.0), color="black"),
    axis.text.x = element_text(size=rel(1.2), color="black")
  )

pdf(file.path(io$outdir,"number_cells_per_omic.pdf"), width=7, height=5)
print(p)
dev.off()
