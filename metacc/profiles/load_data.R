# Sanity checks
# (...)

###########################
## Load methylation data ##
###########################

print("Loading methylation...")

met_list <- list()
# i <- "E8.5_oct20_plate3_A11_L003_CpG"
for (i in opts$met.cells) {

  # met_list[[i]] <- fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,i), select = c(1,2,3,4), colClasses = c("chr"="factor", "start"="integer", "end"="integer", "rate"="numeric")) %>% 
  met_list[[i]] <- fread(sprintf("%s/%s.tsv.gz",io$met_data_raw,i), colClasses = c("chr"="factor", "pos"="integer", "rate"="numeric")) %>% 
    .[,id_met:=i] %>%  # .[,id_met:=factor(i)] %>% 
    .[,c("start","end"):=pos] %>%
    setkey("chr","start","end") %>%
    .[,bp:=start] %>%
    foverlaps(.,anno_df, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    .[,id:=as.character(id)] %>%
    .[,dist:=ifelse(strand %in% c("+","*"),bp-center,center-bp)] %>% 
    .[, dist:=args$met_tile*round(dist/args$met_tile)] %>%
    .[,list(rate=100*mean(rate), N=.N),by=.(id_met,id,dist,anno)]

  print(sprintf("%s: %s features",i,nrow(met_list[[i]])))
}
met.dt <- rbindlist(met_list) %>%
  .[,c("id","context"):=list(as.factor(id),as.factor("CG"))]
  
rm(met_list)

print(object.size(met.dt), units="auto")


#############################
## Load accessibility data ##
#############################

print("Loading accessibility...")

acc_list <- list()
for (i in opts$acc.cells) {

  # acc_list[[i]] <- fread(sprintf("%s/%s.tsv.gz",io$acc_data_raw,i), select = c(1,2,3,4), colClasses = c("chr"="factor", "start"="integer", "end"="integer", "rate"="numeric")) %>% 
  acc_list[[i]] <- fread(sprintf("%s/%s.tsv.gz",io$acc_data_raw,i), colClasses = c("chr"="factor", "pos"="integer", "rate"="numeric")) %>% 
    .[,id_acc:=i] %>% # .[,id_acc:=as.factor(i)] %>% 
    .[,c("start","end"):=pos] %>%
    setkey("chr","start","end") %>%
    .[,bp:=start] %>%
    foverlaps(.,anno_df, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    .[,id:=as.character(id)] %>%
    .[,dist:=ifelse(strand %in% c("+","*"),bp-center,center-bp)] %>% 
    .[, dist:=args$acc_tile*round(dist/args$acc_tile)] %>%
    .[,list(rate=100*mean(rate), N=.N),by=.(id_acc,id,dist,anno)]

  print(sprintf("%s: %s features",i,nrow(acc_list[[i]])))
}
acc.dt <- rbindlist(acc_list) %>%
  .[,c("id","context"):=list(as.factor(id),as.factor("GC"))]
  
rm(acc_list)

print(object.size(acc.dt), units="auto")
