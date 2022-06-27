suppressPackageStartupMessages(library(MOFA2))

###################
## Load settings ##
###################

source(here::here("metacc/mofa/load_settings.R"))
io$outdir <- paste0(io$basedir,"/metacc/mofa")

###############
## Load data ##
###############

source(here::here("metacc/mofa/prepare_data.R"))

#######################
# Create MOFA object ##
#######################

MOFAobject <- create_mofa(data)

# Visualise data structure
plot_data_overview(MOFAobject)

# Data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$use_float32 <- TRUE

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15
model_opts$ard_weights <- FALSE
model_opts$spikeslab_weights <- FALSE

# Training options
train_opts <- get_default_training_options(MOFAobject)

# Prepare MOFA object
MOFAobject <- prepare_mofa(MOFAobject,
  model_options = model_opts,
  data_options = data_opts,
  training_options = train_opts
)

# Train the model
MOFAobject <- run_mofa(MOFAobject, use_basilisk = FALSE)

#########################
## Downstream analysis ##
#########################

# Add sample metadata 
cells <- as.character(unname(unlist(MOFA2::samples_names(MOFAobject))))
samples_metadata(MOFAobject) <- cell_metadata.dt %>% setkey(cell) %>% .[cells]

# Save
saveRDS(MOFAobject, io$outfile)
