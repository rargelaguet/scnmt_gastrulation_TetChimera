create_dir <- function(dir){
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

fwrite_tsv <- partial(fwrite, sep = "\t", na = "NA")


str_from_filename <- function(string, filename){
    split <- strsplit(filename, "_") %>%
        unlist()
    split[grep(toupper(string), toupper(split))][1]
}

rename_stages <- function(filename){
    if (grepl("E4_5-5_5", filename)) {
        return(gsub("E4_5-5_5", "E4.5-5.5", filename))
    }
    
    stage <- strsplit(filename, "_") %>%
        unlist() %>%
        .[grep("E[0-9]$", .)]  %>%
        paste0(".5")
    gsub("E[0-9]_5", stage, filename)
}



# finds bismark files, extracts sample names and creates tsv file names
# does not read or write files
# returns list(infiles, outfiles)



prepare_files <- function(indir, 
                          cg_outdir,
                          gc_outdir = cg_outdir,
                          cg_pattern = "*CpG.cov.gz", 
                          gc_pattern = "*GpC.cov.gz", 
                          exclude_regex = NULL, # exclude files containing this pattern
                          name_index = NULL,
                          name_strings = c("^E[0-9]\\.[0-9]", "Plate|neg", "^[ABCDEFGH][0-9][0-9]$|^[ABCDEFGH][0-9]$"),
                          name_prefix = NULL, # add this to sample names
                          rename_fun = rename_stages
                          ){
    
  in_files <- list(CpG = cg_pattern, GpC = gc_pattern) %>%
        map(~dir(indir, pattern = ., full = TRUE))
    
    if (!is.null(exclude_regex)) {
      in_files <- map(in_files, ~.[!grepl(exclude_regex, .)])
    }
    
    if (is.function(rename_fun)){
        in_files_mod <- map(in_files, map_chr, rename_fun)
    } else {
      in_files_mod <- in_files
    }
    
    out_files <- map(in_files_mod, ~{
        
        if (!is.null(name_index)){
            strsplit(basename(.), "_") %>%
                map(~.[name_index]) %>%
                map_chr(paste, collapse = "_")
        } else {
            map(., function(fn) map_chr(name_strings, str_from_filename, fn)) %>%
                map_chr(paste, collapse = "_")
        }
    })  %>%
        map2(names(.), ~paste0(.x, "_", .y, ".tsv")) %>%
        map2(list(cg_outdir, gc_outdir), ~paste(.y, .x, sep = "/"))
    
    if (!is.null(name_prefix)) {
        out_files <- map(out_files, ~paste0(dirname(.), "/", name_prefix, "_", basename(.)))
    } else {
        
    }
        
    
        list(in_files = in_files, out_files = out_files)
}




# With paired end sequencing, each sample produces 2 bismark files. This 
# function reads in both files into memory and merges
# returns a data.table
merge_reads <- function(files, columns = c(1:2, 5:6)){
  size <- file.size(files)
  if (any(size < 100)) return(data.table(chr = NA, pos = NA, met_reads = NA, nonmet_reads = NA))
  
  map(files, fread, select = columns) %>%
        rbindlist() %>%
        setnames(c("chr", "pos", "met_reads", "nonmet_reads")) %>%
        .[, .(met_reads = sum(met_reads), nonmet_reads = sum(nonmet_reads)), .(chr, pos)] 
}

# Calculates methylation rates and removed non-binary
# input = data.table(chr, pos, met_reads, total_reads)
# output = data.table(chr, pos, met_reads, total_reads, rate)
binarise_methylation_rates <- function(met_dt, binarise = TRUE){
    met_dt[, rate := met_reads / (met_reads + nonmet_reads)]
    setkey(met_dt, rate)
    
    impossible <- met_dt[rate > 1 | rate < 0, .N]
    non_binary <- met_dt[rate > 0 & rate < 1, .N] / met_dt[, .N] 
    
        
    
    
    print(paste("There are ", 
                impossible, 
                " sites with methylation rate higher than 100 or lower than 0"))
    
    print(paste0(round(100 * non_binary, 2), 
                "% sites have non-binary methylation rates"))
    
    if (binarise) met_dt <- met_dt[rate != 0.5][, rate := round(rate)]
    met_dt
}


process_sample <- function(sample, 
                           in_files, 
                           out_files, 
                           context = "CpG",
                           columns = c(1:2, 5:6), 
                           binarise = TRUE,
                           sort = FALSE
                           ){
    if (toupper(context) == "CG") context = "CpG"
    if (toupper(context) == "GC") context = "GpC"
    
    
    short_name <- gsub(".tsv$", "", basename(sample))
    print(paste("Processing", short_name))
    
    index <- which(out_files[[context]] %in% sample)
    in_file <- in_files[[context]][index]
    
    dt <- merge_reads(in_file, columns = columns) %>% 
        binarise_methylation_rates(binarise = binarise)
    
    if (sort) setkey(dt, chr, pos)
    
    fwrite(dt, sample, sep = "\t")
    
    data.table(id = short_name, 
               omic = context, 
               coverage = dt[, .N], 
               mean_rate = dt[, mean(rate) %>% round(2)],
               file = sample
               )
}


process_bismark_files <- function(in_files, 
                                  out_files, 
                                  columns = c(1:2, 5:6), 
                                  binarise = TRUE,
                                  parallel = TRUE, 
                                  sort = TRUE,
                                  gzip = TRUE
                                  ){
    
    # find unique sample names
    out_unique <- map(out_files, unique)
    
    if (parallel) {
        plan(multiprocess)
        map_fun <- function(...){
            future_map(...)
        }
    } else {
        map_fun <- function(...){
            map(...)
        }
    }
        
    # process CpG
    CpG <- map_fun(out_unique$CpG, 
                   process_sample, 
                   in_files, 
                   out_files, 
                   context = "CpG",
                   columns = columns,
                   binarise = binarise,
                   sort = sort
               ) %>%
        rbindlist()
    
    # process GpC
    GpC <- map_fun(out_unique$GpC, 
                   process_sample, 
                   in_files, 
                   out_files, 
                   context = "GpC",
                   columns = columns,
                   binarise = binarise,
                   sort = sort
    ) %>%
        rbindlist()
    
    qc <- rbindlist(list(CpG, GpC))
    
    if (gzip) {
        paste("gzip -f", qc[, file]) %>%
            walk(system)
        qc[, file := paste0(file, ".gz")]
    }
    
    qc
}











