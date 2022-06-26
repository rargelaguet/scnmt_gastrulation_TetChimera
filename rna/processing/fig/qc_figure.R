here::i_am("rna/processing/2_QC.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"processed/rna/metadata.txt.gz")
opts$min_nFeature_RNA <- 4000
opts$max_nFeature_RNA <- 10000
opts$mit_percent_RNA <- 10
opts$rib_percent_RNA <- 20
io$outdir <- file.path(io$basedir,"results/rna/qc/fig"); dir.create(io$outdir, showWarnings = F)

###############
## Load data ##
###############

cell_metadata.dt <- fread(io$metadata) %>% 
    .[!grepl("crispr",class)] %>% 
    .[,class:=ifelse(grepl("WT",class),"WT","Tet-TKO")] %>%
    .[,pass_rnaQC:=nFeature_RNA<=opts$max_nFeature_RNA & nFeature_RNA>=opts$min_nFeature_RNA & mit_percent_RNA<opts$mit_percent_RNA & rib_percent_RNA<opts$rib_percent_RNA]

table(cell_metadata.dt$pass_rnaQC)

##################################
## Plot metrics before QC calls ##
##################################

qc_metrics.to.plot <- c("nFeature_RNA","mit_percent_RNA","rib_percent_RNA")

to.plot <- cell_metadata.dt %>%
    .[mit_percent_RNA<=25] %>% # remove outliers for plotting
    melt(id.vars=c("cell","class"), measure.vars=qc_metrics.to.plot)

tmp <- data.table(
    variable = factor(qc_metrics.to.plot, levels=qc_metrics.to.plot),
    value = c(opts$min_nFeature_RNA, opts$mit_percent_RNA, opts$rib_percent_RNA)
)
to.plot[,variable:=factor(variable, levels=qc_metrics.to.plot)]
facet.labels <- c("nFeature_RNA" = "Num. of genes", "mit_percent_RNA" = "Mitochondrial %", "rib_percent_RNA" = "Ribosomal %")

p <- ggplot(to.plot, aes_string(x="value", fill="class", color="class")) +
    geom_histogram(bins=35) +
    geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
    scale_fill_manual(values=opts$class.colors[unique(to.plot$class)]) + 
    scale_color_manual(values=opts$class.colors[unique(to.plot$class)]) +
    facet_wrap(~variable, scales="free", nrow=1, labeller = as_labeller(facet.labels)) +
    theme_classic() +
    theme(
        axis.text =  element_text(size=rel(0.65)),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        # legend.text = element_text(size=rel(0.85))
        strip.background = element_blank(),
        strip.text = element_text(size=rel(1.15)),
    )

pdf(file.path(io$outdir,"rna_qc_metrics_histogram.pdf"), width=6, height=6)
print(p)
dev.off()

#################################
## Plot metrics after QC calls ##
#################################

to.plot <- cell_metadata.dt %>% .[pass_rnaQC==TRUE] %>%
    melt(id.vars=c("cell","sample","class"), measure.vars=c("nFeature_RNA","mit_percent_RNA","rib_percent_RNA"))

facet.labels <- c("nFeature_RNA" = "Num. of genes", "mit_percent_RNA" = "Mitochondrial %", "rib_percent_RNA" = "Ribosomal %")
to.plot[,variable:=factor(variable, levels=names(facet.labels))]

p <- ggplot(to.plot, aes(x=sample, y=value, fill=class)) +
    ggrastr::geom_jitter_rast(size=0.75, alpha=0.65, width=0.1, shape=21, stroke=0.15) +
    geom_boxplot(outlier.shape=NA, alpha=0.85, coef=1) +
    scale_fill_manual(values=opts$class.colors[unique(to.plot$class)]) + 
    facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(facet.labels)) +
    # scale_fill_manual(values=opts$stage.colors) +
    theme_classic() +
    theme(
        legend.position = "none",
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour="black",size=rel(1.15)),
        axis.text.y = element_text(colour="black",size=rel(1)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=90, hjust=1, vjust=0.5),
        axis.title = element_blank()
    )

pdf(file.path(io$outdir,"rna_qc_metrics_boxplot.pdf"), width=6.5, height=6)
print(p)
dev.off()

