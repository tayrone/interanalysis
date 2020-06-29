# This script aims to identify highly methylated regulons among
# the ones identified as master regulators.

options(stringsAsFactors = F)
library(tidyverse)
library(gdata)
library(ggplot2)
library(RTN)

#---- Loading and preprocessing data ----

load("./interanalysis_files/rdata_files/3_dm_genes.RData")

#dmrs$dmr_id <- as.numeric(dmrs$dmr_id)
#gene_coords$dmr_id <- as.numeric(gene_coords$dmr_id)

load("../expression_analysis/rdata_files/network/g34_rtn.RData")

rm(rtni)

regulons <- tna.get(rtna, what = "regulons")

tfs <- rtna@regulatoryElements

#dm_genes <- gene_coords$hgnc_symbol
#dm_genes <- unique(dm_genes)


#---- Table creation for each regulon ----

# Indexing converts all unexistent values to NA, instead of repeating the vector
#regulons <- sapply(regulons, '[', seq(max(lengths(regulons))))
#regulons <- as.data.frame(regulons)

human_genes_count <- length(rtna@phenotype)

cont_table <- data.frame(matrix(0, nrow = 2, ncol = 2),
                         row.names = c("in_regulon", "out_regulon"))

colnames(cont_table) = c("met", "not_met")

tables_creation <- function(x){
  
  unlist(x)
  
  cont_table["in_regulon", "met"] <- sum(dm_genes %in% x)
  
  
  cont_table["in_regulon", "not_met"] <- sum(!is.na(x)) - sum(dm_genes %in% x)
  
  cont_table["out_regulon", "met"] <- sum(!(dm_genes %in% x))
  
  cont_table["out_regulon", "not_met"] <- 
    human_genes_count - sum(cont_table)
  
  return(cont_table)
  
}


# hm_regulon_test <- function(object){
#   #cat(as.name(object), "Expected: \n")
#   #print(chisq.test(get(object))$expected)
#   
#   if(any(chisq.test(object)$expected < 5)){
#     #cat("Fisher test: \n")
#     #print(fisher.test(get(object)))
#     #cat("P-value: \n")
#     return(fisher.test(object)$p.value)
#     
#   }else{
#     #cat("Chi test: \n")
#     #print(chisq.test(get(object)))
#     #cat("P-value: \n")
#     return(chisq.test(object)$p.value)
#   }
# }


#---- Significance test ----

tables <- lapply(regulons, tables_creation)

p_values <- lapply(tables, fisher.test)

p_values <- sapply(p_values, function(x) return(x$p.value))

adjusted_p <- p.adjust(p_values, method = "bonferroni")

hm_regulons <- names(adjusted_p)[adjusted_p < 0.01]

#----

gdata::keep(dm_genes, hm_regulons, adjusted_p, tables, tfs, sure = T)

save.image("./interanalysis_files/rdata_files/5_dm_regulons.RData")
