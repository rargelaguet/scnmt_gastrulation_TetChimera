library(RColorBrewer)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/human_embryo_multiomics/settings.R")
  source("/Users/ricard/human_embryo_multiomics/met/differential/feature_level/lineages/analysis/utils.R")
} else {
  stop("Computer not recognised")
}


#########
## I/O ##
#########

io$indir <- paste0(io$basedir,"/met/results/differential/feature_level/lineages/old")
io$outdir <- paste0(io$basedir,"/met/results/differential/feature_level/lineages/old/pdf")

#############
## Options ##
#############

opts$comparisons <- list(
  "PrE_commitment_vs_TE_commitment" = c("ICM_vs_PrE","ICM_vs_Epiblast")
)

# Select genomic contexts
opts$annos <- NULL
if (is.null(opts$annos)) {
  opts$annos <- list.files(io$met_data_parsed, pattern=".tsv.gz") %>% gsub(".tsv.gz","",.)
}

###############
## Load data ##
###############

# Load precomputed differential results
diff.results <- lapply(names(opts$comparisons), function(i) 
  lapply(opts$annos, function(j) {
    i <- opts$comparisons[[i]]
    file1 <- sprintf("%s/%s_%s.txt.gz",io$indir,i[[1]],j)
    file2 <- sprintf("%s/%s_%s.txt.gz",io$indir,i[[2]],j)
    if (file.exists(file1) & file.exists(file2)) {
      foo <- fread(file1) %>% .[,subcomparison:=i[[1]]]
      bar <- fread(file2)  %>% .[,subcomparison:=i[[2]]]
      foobar <- rbind(foo,bar) %>% .[,anno:=as.factor(j)]
      return(foobar)
    }
  }) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist

diff.results <- diff.results[,c("id","diff","log_padj_fdr","anno","subcomparison","comparison")]

##############
## Barplots ##
##############

to.plot <- diff.results %>% split(.$comparison) %>% map(function(x)
  x %>% dcast(id+anno~subcomparison, value.var="diff") %>% .[complete.cases(.)] 
) %>% rbindlist %>% setnames(c("id","anno","diff_comparison1","diff_comparison2"))

p <- ggplot(to.plot, aes(x=diff_comparison1, y=diff_comparison2)) +
  facet_wrap(~anno, scales="fixed") +
  # geom_point() +
  geom_bin2d(bins=100, fill="black") +
  geom_smooth(method="lm") +
  scale_fill_viridis_c() +
  facet_wrap(~anno) +
  theme_pubr() +
  labs(x="Morula vs ICM (Methylation %)", y="Morula vs TE_mural (Methylation %)") +
  theme(
    axis.text = element_text(size=rel(0.75))
  )

# pdf(sprintf("%s/scatterplot_%s_AND_%s.pdf",io$outdir,opts$comparisons[[1]][1], opts$comparisons[[1]][2]), width = 13, height = 13, useDingbats = FALSE)
print(p)
# dev.off()