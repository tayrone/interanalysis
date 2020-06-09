options(stringsAsFactors = F)
library(tidyverse)
library(gdata)

load("./interanalysis_files/rdata_files/3_dmgs_to_genes.RData")

dmrs$dmr_id <- as.numeric(dmrs$dmr_id)
gene_coords$dmr_id <- as.numeric(gene_coords$dmr_id)

load("./expression_analysis/rdata_files/network_g34_rtn.RData")

rm(rtni)

regulons <- tna.get(rtna, what = "regulons")

tfs <- rtna@regulatoryElements

#----

# Indexing converts all unexistent values to NA, instead of repeating the vector
#regulons <- sapply(regulons, '[', seq(max(lengths(regulons))))
#regulons <- as.data.frame(regulons)

dm_genes <- gene_coords$hgnc_symbol
dm_genes <- unique(dm_genes)


cont_table <- data.frame(matrix(0, nrow = 2, ncol = 2),
                         row.names = c("in_regulon", "out_regulon"))

colnames(cont_table) = c("dm", "not_dm")


tables_creation <- function(x){
  
  unlist(x)
  
  cont_table["in_regulon", "dm"] <- sum(dm_genes %in% x)
  cont_table["out_regulon", "dm"] <- sum(!(dm_genes %in% x))
  
  cont_table["in_regulon", "not_dm"] <- sum(!is.na(x)) - sum(dm_genes %in% x)
  
  cont_table["out_regulon", "not_dm"] <- 
    length(unique(c(tfs, dm_genes))) - sum(cont_table)
  
  return(cont_table)
  
}

tables <- lapply(regulons, tables_creation)


sig_test <- function(x){
  
  x <- as.data.frame(x)
  
  #cat(as.name(object), "Expected: \n")
  
  #print(chisq.test(get(object))$expected)
  
  if(any(chisq.test(x)$expected < 5)){
    
    #cat("Fisher test: \n")
    #print(fisher.test(get(object)))
    
    #cat("P-value: \n")
    return(fisher.test(x)$p.value)
    
  }else{
    
    #cat("Chi test: \n")
    #print(chisq.test(get(object)))
    
    #cat("P-value: \n")
    return(chisq.test(x)$p.value)
    
  }
  #cat("\n----------\n")
}

p_values <- sapply(tables, sig_test)


length(unique(c(tfs, dm_genes)))

names(p_values[p_values <= 0.1])


#----


in_and_dm <- sapply(tables, function(x){ 
                              x <- as.data.frame(x)
                              x["in_regulon", "dm"]})

above_avg <- function(x){ 
  x <- as.data.frame(x)
  if(x["in_regulon", "dm"] > mean(in_and_dm)){
    return(T) 
  }else{
    return(F)
  } 
}


dm_regulons <- sapply(tables, above_avg)

dm_regulons <- names(tables)[dm_regulons]


gdata::keep(dmrs, gene_coords, dm_regulons, p_values, sure = T)

save.image("./interanalysis_files/rdata_files/dm_regulons.RData")
