suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(MOFA2))

############################
## Define I/O and options ##
############################

source("/Users/ricard/scnmt_gastrulation_TetChimera/acc/dimensionality_reduction/load_settings.R")

################
## Load model ##
################

# file <- paste0(io$basedir,"/mofa/hdf5/model_prom_2000_2000_E5_E6.hdf5")
# model <- load_model(file)

##################
## Rename views ##
##################

# opts$rename.views <- c(
#   "acc_Promoters" = "Promoter accessibility",
#   # "acc_Gene bodies" = "Genebody accessibility",
#   # "acc_E3.5 enhancers" = "E3.5 enhancer accessibility",
#   "acc_E7.5 enhancers" = "Enhancer accessibility",
#   # "acc_E10.5 enhancers" = "E10.5 enhancer accessibility",
#   
#   "RNA$" = "RNA expression"
# )
# 
# views(model) = stringr::str_replace_all(views(model), opts$rename.views)
# 
# model <- subset_views(model, views=unname(opts$rename.views) )

#########################
## Add sample metadata ##
#########################

cells <- as.character(unname(unlist(MOFA2::samples_names(model))))

sample_metadata_filt <- sample_metadata %>% setkey(id_acc) %>% .[cells] %>%
  .[,sample:=id_acc]

stopifnot(all(cells==sample_metadata_filt$sample))
samples_metadata(model) <- sample_metadata_filt

####################
## Subset factors ##
####################

# r2 <- model@cache$variance_explained$r2_per_factor
# factors <- sapply(r2, function(x) x[,"RNA"]>0.01)
# model <- subset_factors(model, which(apply(factors,1,sum)>=1))
# factors(model) <- paste("Factor",1:get_dimensions(model)[["K"]], sep=" ")

acc <- colMeans(model@data[[1]][[1]],na.rm=T)
acc <- 100*2**acc/(1+2**acc)

plot_variance_explained(model)

plot_factors(model, factors = 1:3, color_by = "stage")

plot_factor(model, factors = 1, color_by = acc)

plot_weights(model, view=1, factor = 1, nfeatures = 0, scale = F)

# plot_data_scatter(model, factor=1, view=1, color_by = "lab", features = 8, dot_size = 2)


