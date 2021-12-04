###########################
## Load default settings ##
###########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation_TetChimera//settings.R")
} else {
  stop("Computer not recognised")
}

io$indir <- paste0(io$basedir,"/met/differential/feature_level/pseudobulk")
io$outdir <- paste0(io$basedir,"/met/differential/feature_level/pseudobulk/pdf")

#############
## Options ##
#############

opts$comparisons <- c(
  # "ICM_vs_TE_mural"
  "Epiblast_vs_PrE"
)

# Select genomic contexts
opts$annos <- NULL
if (is.null(opts$annos)) {
  opts$annos <- list.files(io$features.dir, pattern=".bed.gz") %>% gsub(".bed.gz","",.)
}

###############
## Load data ##
###############


# Load precomputed differential results
diff.results <- lapply(opts$comparisons, function(i) 
  lapply(opts$annos, function(j) {
    file <- sprintf("%s/%s_%s.txt.gz",io$indir,i,j)
    if (file.exists(file)) {
      fread(file) %>% .[,anno:=as.factor(j)]
    } else {
      cat(sprintf("%s does not exist",file))
    }
  }
  ) %>% rbindlist %>% .[,comparison:=i] 
) %>% rbindlist# %>% .[complete.cases(.)]

diff.results[,diff:=-diff]

##############
## Barplots ##
##############

gg_barplot <- function(tmp, title = "", ylim=NULL) {
  
  if (is.null(ylim)) {
    ylim <- c(min(tmp$value, na.rm=T), max(tmp$value, na.rm=T))
  }
  
  p <- ggplot(tmp, aes(x=anno, y=value, group=anno)) +
    geom_bar(aes(fill=anno), color="black", stat="identity", position="dodge") +
    geom_hline(yintercept=0, color="black") +
    scale_y_continuous(limits=c(ylim[1],ylim[2])) +
    labs(title="", x="", y="Fraction of hits") +
    theme_classic() +
    theme(
      plot.title = element_text(size=11, face='bold', hjust=0.5),
      axis.text = element_text(size=rel(1.0), color='black'),
      # axis.text.x = element_blank(),
      # axis.text.x = element_text(colour="black", angle=30, vjust=1, hjust=1),
      axis.text.x = element_text(colour="black",size=rel(0.8), angle=90, hjust=1, vjust=0.5),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size=rel(1.0), color='black'),
      axis.line = element_line(color="black"),
      legend.position="none"
    )
  
  return(p)
}

to.plot <- diff.results %>%
  .[,.(number_positive_hits=mean(diff_rate>25), number_negative_hits=mean(diff_rate<(-25))), by=c("anno","comparison")] %>%
  .[,number_negative_hits:=-number_negative_hits] %>%
  melt(id.vars=c("anno","comparison"))

ylim <- c(min(to.plot$value), max(to.plot$value))

for (i in unique(diff.results$comparison)) {
  p <- gg_barplot(to.plot[comparison==i], title=i, ylim=ylim)
  
  pdf(sprintf("%s/barplots_%s.pdf",io$outdir,i), width=12, height=4)
  print(p)
  dev.off()
}
