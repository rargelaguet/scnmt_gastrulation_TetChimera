# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$trajectory.dir <- file.path(io$basedir,"results_new/rna/trajectories/blood")
io$outdir <- file.path(io$basedir,"results_new/rna/trajectories/blood/pdf")
dir.create(io$outdir, showWarnings=F)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(file.path(io$trajectory.dir,"blood_sample_metadata.txt.gz")) 

#####################
## Load trajectory ##
#####################

trajectory.dt <- fread(file.path(io$trajectory.dir,"blood_trajectory.txt.gz")) %>%
  setnames(c("DC1","DC2"),c("V1","V2"))

##################
## Load Met/Acc ##
##################

global_rates_metacc.dt <- sample_metadata.dt[,c("cell","id_rna","met_rate","acc_rate")]# %>% 
  # setnames(c("cell","CG","GC"))# %>%
  # melt(id.vars="cell", variable.name="context", value.name="global_rate")

# met.dt <- fread(args$met_file) %>% 
#   setnames(c("id_met","id","anno","Nmet","Ntotal","rate")) %>% 
#   .[id_met%in%opts$met_cells]
# 
# acc.dt <- fread(args$acc_file) %>% 
#   setnames(c("id_acc","id","anno","Nmet","Ntotal","rate")) %>% 
#   .[id_acc%in%opts$acc_cells]
# 
# # Add common cell identifier
# met.dt <- merge(met.dt, sample_metadata[,c("cell","id_met")], by="id_met") %>% 
#   .[,id_met:=NULL] %>% .[,context:="CG"]
# acc.dt <- merge(acc.dt, sample_metadata[,c("cell","id_acc")], by="id_acc") %>% 
#   .[,id_acc:=NULL] %>% .[,context:="GC"]
# 
# # Merge
# metacc.dt <- rbind(met.dt,acc.dt)

#####################################
## Impute missing methylation data ##
#####################################

mtx <- global_rates_metacc.dt %>% .[,c("id_rna","met_rate")] %>% matrix.please %>% t
trajectory.mtx <- trajectory.dt[,c("id_rna","V1","V2")] %>% matrix.please
stopifnot(colnames(mtx)==rownames(trajectory.mtx))
global_rates_metacc.dt[,met_rate_imputed:=smoother_aggregate_nearest_nb(mat=mtx, D=pdist(trajectory.mtx), k=25)]
  
#########################
## Load RNA expression ##
#########################

sce <- load_SingleCellExperiment(io$rna.sce, normalise = T, cells=sample_metadata.dt$id_rna)

# Add sample metadata as colData
colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("id_rna") %>% DataFrame

###################################
## Plot dimensionality reduction ##
###################################

to.plot <- trajectory.dt %>% merge(sample_metadata.dt, by="id_rna")

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=celltype)) +
  # geom_point(size=3, shape=21, stroke=0.25) +
  geom_jitter(size=3, shape=21, stroke=0.25, width=0.002, height=0.002) +
  scale_fill_manual(values=opts$celltype.colors) +
  ggrepel::geom_text_repel(aes_string(label="celltype"), data=to.plot[,.(V1=median(V1), V2=median(V2)), by="celltype"]) +
  theme_classic() +
  labs(x="", y="") +
  ggplot_theme_NoAxes() +
  theme(
    legend.position="none"
  )

pdf(paste0(io$outdir,"/trajectory_coloured_by_cellype.pdf"), width=6, height=5)
print(p)
dev.off()

##########################
## Plot gene expression ##
##########################

genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3")

for (i in genes.to.plot) {
  
  to.plot <- data.table(
    id_rna = colnames(sce),
    expr = logcounts(sce)[i,]
  ) %>% merge(trajectory.dt, by="id_rna")
  
  p <- ggplot(to.plot, aes(x=V1, y=V2, fill=expr)) +
    # geom_point(size=3, shape=21, stroke=0.2) +
    geom_jitter(size=3, shape=21, stroke=0.25, width=0.002, height=0.002) +
    scale_fill_gradient2(low = "gray50", mid="gray90", high = "darkgreen") +
    # ggrepel::geom_text_repel(aes_string(label="celltype"), data=to.plot[,.(V1=median(V1), V2=median(V2)), by="celltype"]) +
    theme_classic() +
    labs(x="", y="", title=i) +
    ggplot_theme_NoAxes() +
    theme(
      legend.position="none"
    )
  
  pdf(file.path(io$outdir,sprintf("trajectory_coloured_by_%s_expr.pdf",i)), width=6, height=5)
  print(p)
  dev.off()
    
}


##################
## Plot Met/Acc ##
##################
  
to.plot <- global_rates_metacc.dt %>% merge(trajectory.dt, by="id_rna") %>% 
  merge(sample_metadata.dt[,c("id_rna","class")],by="id_rna")

# Just for visualisation purposes
max.meth <- 83
to.plot[met_rate_imputed>max.meth,met_rate_imputed:=max.meth]

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=met_rate_imputed)) +
  geom_jitter(size=3, shape=21, stroke=0.25, width=0.002, height=0.002) +
  theme_classic() +
  scale_fill_distiller(palette = "YlOrRd", direction=1) +
  # labs(x="", y="", title="Global DNA methylation") +
  ggplot_theme_NoAxes() +
  theme(
    legend.position="none"
  )

pdf(file.path(io$outdir,"trajectory_coloured_by_global_met.pdf"), width=7, height=6)
print(p)
dev.off()
  
