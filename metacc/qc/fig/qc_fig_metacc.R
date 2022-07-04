here::here("metacc/qc/qc.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$bisulfite_conversion <- file.path(io$basedir,"original/conversion.tsv.gz") 

# Options
opts$context <- "CG"

if (opts$context=="CG") {
  io$metadata <- file.path(io$basedir,"results/met/stats/sample_metadata_after_met_stats.txt.gz")
  io$outdir <- file.path(io$basedir,"results/met/qc/fig"); dir.create(io$outdir, showWarnings = F)
  opts$minimum_number_sites <- 5000
  opts$min_rate <- 50
  opts$max_rate <- 100
} else if (opts$context=="GC") {
  io$metadata <- file.path(io$basedir,"results/acc/stats/sample_metadata_after_acc_stats.txt.gz")
  io$outdir <- file.path(io$basedir,"results/acc/qc/fig"); dir.create(io$outdir, showWarnings = F)
  opts$minimum_number_sites <- 10000
  opts$min_rate <- 10
  opts$max_rate <- 40
}

###################
## Load metadata ##
###################

cell_metadata.dt <- fread(io$metadata) %>%
  .[,class:=ifelse(grepl("WT",class),"WT","Tet-TKO")]

##########################
## Plot before QC calls ##
##########################


##############
## QC calls ##
##############

if (opts$context=="CG") {
  cell_metadata.dt %>% .[,pass_metQC:=met_rate<=opts$max_rate & met_rate>=opts$min_rate & nCG>=opts$minimum_number_sites]
} else if (opts$context=="GC") {
  cell_metadata.dt %>% .[,pass_accQC:=acc_rate<=opts$max_rate & acc_rate>=opts$min_rate & nGC>=opts$minimum_number_sites]
}

#############################
## Boxplots after QC calls ##
#############################

if (opts$context=="CG") {
  to.plot <- cell_metadata.dt %>% .[pass_metQC==TRUE] %>% setnames(c("nCG","met_rate"),c("N","rate"))
} else if (opts$context=="GC") {
  to.plot <- cell_metadata.dt %>% .[pass_accQC==TRUE] %>% setnames(c("nGC","acc_rate"),c("N","rate"))
}

to.plot.melted <- to.plot %>% 
  .[,log10_N:=log10(N)] %>%
  melt(id.vars=c("cell","sample","class"), measure.vars=c("log10_N","rate")) 

p <- ggplot(to.plot.melted, aes(x=sample, y=value, fill=class)) +
  ggrastr::geom_jitter_rast(size=1, alpha=0.65, width=0.1, shape=21, stroke=0.15) +
  geom_boxplot(outlier.shape=NA, alpha=0.85, coef=1) +
  scale_fill_manual(values=opts$class.colors[unique(to.plot$class)]) + 
  facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(c("log10_N"="Num. of observations", "rate"="Rate (%)"))) +
  # scale_fill_manual(values=opts$stage.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(colour="black",size=rel(1.25)),
    axis.text.y = element_text(colour="black",size=rel(1.25)),
    axis.text.x = element_text(colour="black",size=rel(0.75), angle=90, hjust=1, vjust=0.5),
    axis.title = element_blank()
  )

pdf(file.path(io$outdir,sprintf("%s_qc_metrics_boxplot.pdf",opts$context)), width=7.5, height=7)
print(p)
dev.off()

################################
## Scatterplot after QC calls ##
################################

if (opts$context=="CG") {
  to.plot <- cell_metadata.dt %>% .[!is.na(id_met)] %>% 
    setnames(c("nCG","met_rate"),c("N","rate")) %>% setnames("pass_metQC","pass_QC")
} else if (opts$context=="GC") {
  to.plot <- cell_metadata.dt %>% .[!is.na(id_acc)] %>% 
    setnames(c("nGC","acc_rate"),c("N","rate")) %>% setnames("pass_accQC","pass_QC")
}

p <- ggplot(to.plot, aes(x=rate, y=log10(N), fill=pass_QC)) +
  geom_point(shape=21, size=2.5, stroke=0.2) +
  geom_hline(yintercept=log10(opts$minimum_number_sites), linetype="dashed", color="black") +
  geom_vline(xintercept=opts$min_rate, linetype="dashed", color="black") +
  geom_vline(xintercept=opts$max_rate, linetype="dashed", color="black") +
  labs(x="Rate (%)", y="Number of observations") +
  scale_fill_manual(values=c("FALSE"="gray70", "TRUE"="red")) +
  theme_classic() +
  theme(
    axis.text = element_text(color="black"),
    legend.position = "none"
  )

pdf(file.path(io$outdir,sprintf("%s_qc_metrics_scatterplot.pdf",opts$context)), width=6, height=6)
print(p)
dev.off()


################################
## Bisulfite conversion rates ##
################################

bisulfite_conversion_rates.dt <- fread(io$bisulfite_conversion) %>%
  .[,cell:=gsub("_L00.","",cell)] %>% 
  merge(cell_metadata.dt[pass_metQC==T,c("cell","method","stage","class","sample")])

unique(bisulfite_conversion_rates.dt$cell)

to.plot <- bisulfite_conversion_rates.dt[,.(conversion=mean(conversion), sd=sd(conversion)), by=c("class","sample")]

p <- ggplot(to.plot, aes(x=sample, y=conversion, fill=class)) +
  # ggrastr::geom_jitter_rast(size=1, alpha=0.65, width=0.1, shape=21, stroke=0.15) +
  geom_bar(position=position_dodge(), stat="identity", color="black") +
  geom_errorbar(aes(ymin=conversion-sd, ymax=conversion+sd), width=0.4, colour="black", alpha=0.9, size=0.5) +
  scale_fill_manual(values=opts$class.colors[unique(bisulfite_conversion_rates.dt$class)]) + 
  labs(x="", y="Bisulfite conversion rate (%)") +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(1)),
    axis.text.x = element_text(colour="black",size=rel(0.75), angle=90, hjust=1, vjust=0.5),
    axis.title = element_blank()
  )

pdf(file.path(io$outdir,"bisulfite_conversion_rates.pdf"), width=5.5, height=6)
print(p)
dev.off()