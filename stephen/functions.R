read_covs <- function(file1, file2 = NULL, name = "", cols = c("chr", "pos", "met_reads", "nonmet_reads", "rate")){
  dt <- map(c(file1, file2), ~fread(cmd = paste("zcat <", .), select = c(1:2, 5:6))) %>%
    rbindlist() %>%
      .[, .(met_reads = sum(V5), nonmet_reads = sum(V6)), .(V1, V2)] %>%
        .[, .(chr = V1,
              pos = V2,
              met_reads,
              nonmet_reads,
              rate = met_reads / (met_reads + nonmet_reads))] 
  
  dt[, .SD, .SDcol = cols]
}

annotate_met <- function(raw_met, anno){
  if ("pos" %in% colnames(raw_met)) raw_met[, c("start", "end") := pos]
  
  global_rate <- raw_met[, mean(rate)]
  
  setkey(anno, chr, start, end)
  setkey(raw_met, chr, start, end)
  
  foverlaps(raw_met, anno, nomatch = 0L) %>%
    .[, .(rate = mean(rate), .N, global_rate = global_rate), .(id, anno)]
}

fread_gz <- function(file, ...){
  fread(cmd = paste("zcat <", file), ...)
}

fwrite_tsv <- partial(fwrite, sep = "\t", na = "NA")

sce_to_dt <- function(sce){
  
  dt=logcounts(sce) %>%
    as.data.frame() %>%
      setDT(keep.rownames = "ens_id") %>%
        melt(id.vars = "ens_id", variable.name = "id_rna", value.name = "expr")
  
  
}


pseudobulk_met_old <- function(files){
  map(files, fread_gz, select = c(1:2, 5)) %>%
    rbindlist() %>%
      setkey(chr, pos) %>%
       .[, .(rate = mean(rate), .N), .(chr, pos)]
}
fread_and_sum <- function(dt, file){
  dt2 <- fread_gz(file, select = c(1:2, 5)) 
  dt2[, c("met_reads", "nonmet_reads", "rate") := .(rate, 1-rate, NULL)]
  rbind(dt, dt2) %>%
    .[, .(met_reads = sum(met_reads),
          nonmet_reads = sum(nonmet_reads)),
      .(chr, pos)]
}
pseudobulk_met<- function(files){
  init <- data.table(NULL)
  dt <- purrr::reduce(files, fread_and_sum, .init = init) 
  dt[, rate := met_reads / (met_reads + nonmet_reads)]
}



find_nearby_genes <- function(anno, gene_meta, dist_cutoff = 1e5){
  if ("symbol" %in% colnames(gene_meta)) setnames(gene_meta, "symbol", "gene")
  gene_meta[, chr := gsub("chr", "", chr)]
  gene_meta[strand == "+", tss := start]
  gene_meta[strand == "-", tss := end]
  gene_meta[, c("start", "end") := .(start - dist_cutoff, end + dist_cutoff)]
  
  setkey(gene_meta, chr, start, end)
  setkey(anno, chr, start, end)
  
  foverlaps(anno, gene_meta) %>%
    .[, tss_dist := tss - i.start + (i.end-i.start)/2] %>%
      .[, .(id, anno, ens_id, gene, tss_dist)]
  
}

wtd_cor <- function(x, y, wtd){
  na_return <- list(cor = as.numeric(NA), p = as.numeric(NA))
  if (length(!is.na(x)) < 10 || length(!is.na(y)) < 10) return(na_return)
  if (var(x, na.rm = TRUE) == 0 || var(y, na.rm = TRUE) == 0) return(na_return)
  w <- weights::wtd.cor(x, y, wtd)
  list(cor = w[1], p = w[4])
}


correlate_met_rna <- function(met_dt, 
                              rna_dt, 
                              anno, 
                              gene_meta, 
                              min_cells = 20,
                              var_cutoff = 0.5, # keep this fraction (i.e. 0.25 = top 25%)
                              dist_cutoff = 1e5){
  
  genes <- find_nearby_genes(anno, gene_meta, dist_cutoff = dist_cutoff)
  rna_dt[, ncells := .N, gene]
  met_dt[, ncells := .N, .(id, anno)]
  
  rna_keep <- rna_dt[ncells > min_cells, var(exp), gene] %>%
    .[order(-rank(V1))] %>%
      .[1: (.N * var_cutoff), gene]
  
  met_keep <- met_dt[ncells > min_cells, var(rate), .(id)] %>%
    .[order(-rank(V1))] %>%
      .[1: (.N * var_cutoff), id]
  
  cors <- merge(met_dt[id %in% met_keep], genes, by = c("anno", "id"), allow.cartesian = TRUE) %>%
    merge(rna_dt[gene %in% rna_keep], by = c("gene", "sample"), allow.cartesian = TRUE) %>%
      .[, wtd_cor(exp, rate, N), .(gene, ens_id, anno, id, tss_dist)]
  
  cors[, padj := p.adjust(p, method = "fdr")]
  cors[, logpadj := -log10(padj)]
  cors[order(padj)]
}

