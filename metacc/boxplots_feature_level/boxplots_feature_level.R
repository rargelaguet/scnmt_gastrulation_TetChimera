here::here("metacc/boxplots_feature_level/boxplots_feature_level.R")

###################
## Load settings ##
###################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$outdir <- file.path(io$basedir,"results/metacc/boxplots_feature_level")

# Options
opts$met.annos <- c(
  "prom_2000_2000" = "Promoters"
  )

opts$acc.annos <- c(
  "prom_200_200" = "Promoters"
)


###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata)

opts$met_cells <- sample_metadata[pass_metQC==TRUE,id_met]
opts$acc_cells <- sample_metadata[pass_accQC==TRUE,id_met]

###############
## Load data ##
###############

met_dt <- lapply(names(opts$met.annos), function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,n)) %>% .[V1%in%opts$met_cells]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","Ntotal","rate"))

acc_dt <- lapply(names(opts$acc.annos), function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$acc_data_parsed,n)) %>% .[V1%in%opts$acc_cells]
}) %>% rbindlist %>% setnames(c("id_acc","id","anno","Nmet","Ntotal","rate"))

###########
## Merge ##
###########

# Add common cell ID
acc_dt <- merge(acc_dt, sample_metadata[,c("cell","id_acc")], by="id_acc") %>% .[,id_acc:=NULL]
met_dt <- merge(met_dt, sample_metadata[,c("cell","id_met")], by="id_met") %>% .[,id_met:=NULL]

# Add genomic context ID
# met_dt[,anno:=factor("DNA methylation")]
# acc_dt[,anno:=factor("Chr. accessibility")]

# Rename annotations
met_dt <- met_dt %>% .[,anno:=stringr::str_replace_all(anno,opts$met.annos)]
acc_dt <- acc_dt %>% .[,anno:=stringr::str_replace_all(anno,opts$acc.annos)]

# Merge
to.plot <- merge(met_dt,acc_dt, by=c("cell","id","anno"), suffixes = c("_met", "_acc")) %>% 
  melt(id.vars=c("cell","id","anno")) %>%
  merge(sample_metadata[,c("cell","sample","celltype.mapped_mnn")], by="cell") 

##########################
## Plot absolute levels ##
##########################

p <- ggplot(to.plot, aes(x = sample, y = rate, fill = class)) +
  facet_wrap(~anno)
  geom_boxplot(outlier.shape=NA, coef=1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.1, alpha = 0.2) +
  labs(x="", y="Methylation") +
  facet_wrap(~anno, scales="fixed") +
  #scale_fill_manual(values=opts$colors) +
  # coord_cartesian(ylim=c(5,60)) +
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  theme(
    #legend.position = "none",
    axis.text = element_text(color="black")
  )
p


##########################
## Plot relative levels ##
##########################

# Normalise rates by genome-wide levels
# sample_metadata[,c("")]
