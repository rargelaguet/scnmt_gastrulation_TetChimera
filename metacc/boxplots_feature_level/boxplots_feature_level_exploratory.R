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
opts$acc_cells <- sample_metadata[pass_accQC==TRUE,id_acc]

###############
## Load data ##
###############

met.dt <- lapply(names(opts$met.annos), function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,n)) %>% .[V1%in%opts$met_cells]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","Ntotal","rate"))

acc.dt <- lapply(names(opts$acc.annos), function(n) {
  fread(sprintf("%s/%s.tsv.gz",io$acc_data_parsed,n)) %>% .[V1%in%opts$acc_cells]
}) %>% rbindlist %>% setnames(c("id_acc","id","anno","Nmet","Ntotal","rate"))

###########
## Merge ##
###########

# Add common cell identifier
met.dt <- merge(met.dt, sample_metadata[,c("cell","id_met")], by="id_met") %>% 
  .[,id_met:=NULL] %>% .[,context:="CG"]
acc.dt <- merge(acc.dt, sample_metadata[,c("cell","id_acc")], by="id_acc") %>% 
  .[,id_acc:=NULL] %>% .[,context:="GC"]

# Rename annotations
met.dt <- met.dt %>% .[,anno:=stringr::str_replace_all(anno,opts$met.annos)]
acc.dt <- acc.dt %>% .[,anno:=stringr::str_replace_all(anno,opts$acc.annos)]

# Merge
metacc.dt <- rbind(met.dt,acc.dt)

###############################
## Prepare data for plotting ##
###############################

global_rates.dt <- sample_metadata[,c("cell","met_rate","acc_rate")] %>% 
  setnames(c("cell","CG","GC")) %>%
  melt(id.vars="cell", variable.name="context", value.name="global_rate")

# Normalise rates by genome-wide levels
to.plot <- metacc.dt %>%
  .[,.(Nmet=sum(Nmet), Ntotal=sum(Ntotal)), by=c("cell","anno","context")] %>% 
  .[Ntotal>=250] %>%
  .[,rate:=100*(Nmet/Ntotal)] %>%
  merge(sample_metadata[,c("cell","sample","celltype.mapped")], by="cell") %>%
  merge(global_rates.dt, by=c("cell","context")) %>%
  .[,rate_norm:=rate/global_rate]

##########################
## Plot absolute levels ##
##########################

for (i in unique(to.plot$anno)) {
  p <- ggplot(to.plot[anno==i], aes(x = sample, y = rate, fill=context)) +
    # facet_wrap(~context) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.5, alpha = 0.4, shape=21) +
    facet_wrap(~context, scales="free_y") +
    labs(x="", y="Met/Acc levels (%)") +
    scale_color_manual(values=opts$context.colors) +
    scale_fill_manual(values=opts$context.colors) +
    # coord_cartesian(ylim=c(5,60)) +
    theme_bw() +
    guides(x = guide_axis(angle = 90)) +
    theme(
      legend.position = "none",
      axis.text = element_text(color="black")
    )
  
  pdf(file.path(io$outdir,sprintf("boxplots_metacc_%s.pdf",i)), width=8, height=6)
  print(p)
  dev.off()
}


##########################
## Plot relative levels ##
##########################

for (i in unique(to.plot$anno)) {
  p <- ggplot(to.plot[anno==i], aes(x = sample, y = rate_norm, fill=context)) +
    # facet_wrap(~context) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 0.5, alpha = 0.4, shape=21) +
    facet_wrap(~context, scales="free_y") +
    geom_hline(yintercept = 1, linetype="dashed") +
    labs(x="", y="Met/Acc levels (relative to background)") +
    scale_color_manual(values=opts$context.colors) +
    scale_fill_manual(values=opts$context.colors) +
    # coord_cartesian(ylim=c(5,60)) +
    theme_bw() +
    guides(x = guide_axis(angle = 90)) +
    theme(
      legend.position = "none",
      axis.text = element_text(color="black")
    )
  
  pdf(file.path(io$outdir,sprintf("boxplots_metacc_relative_%s.pdf",i)), width=8, height=6)
  print(p)
  dev.off()
}
