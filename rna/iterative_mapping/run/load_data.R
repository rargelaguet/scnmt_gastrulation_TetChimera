
################
## Load atlas ##
################

# Load cell metadata
meta_atlas <- fread(io$atlas.metadata) %>%
  .[stripped==F & doublet==F & stage%in%args$atlas_stages] %>%
  .[,celltype:=stringr::str_replace_all(celltype,"_", " ")]

if (isTRUE(args$test)) meta_atlas <- head(meta_atlas,n=1000)

# Load SingleCellExperiment
sce_atlas  <- readRDS(io$atlas.sce)[,meta_atlas$cell]
# sce_atlas$celltype <- sce_atlas$celltype 

# Add metadata to the SCE object
colData(sce_atlas) <- meta_atlas %>% tibble::column_to_rownames("cell") %>% DataFrame

################
## Load query ##
################

# Load cell metadata
meta_query <- fread(io$metadata) %>% 
  .[pass_rnaQC==T & plate%in%args$query_plates] %>%
  setnames("id_rna","cell")

# Load SingleCellExperiment
sce_query <- readRDS(io$sce)[,meta_query$cell]

#############
## Prepare ## 
#############

# Filter out non-expressed genes
sce_query <- sce_query[rowMeans(counts(sce_query))>1e-5,]
sce_atlas <- sce_atlas[rowMeans(counts(sce_atlas))>1e-5,]

# Intersect genes
genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

# Load gene markers to be used as HVGs
marker_genes.dt <- fread(io$atlas.marker_genes)
marker_genes.dt <- marker_genes.dt[,head(.SD,n=50),by="celltype"]
marker_genes <- unique(marker_genes.dt$ens_id)
marker_genes <- marker_genes[marker_genes%in%genes.intersect]

stopifnot(all(marker_genes%in%rownames(sce_atlas)))
stopifnot(all(marker_genes%in%rownames(sce_query)))

# Update SingleCellExperiment objects
sce_query <- sce_query[marker_genes,]
sce_atlas <- sce_atlas[marker_genes,]
