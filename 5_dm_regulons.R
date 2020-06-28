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

#---- Significance test ----

human_genes_count <- length(rtna@phenotype)
dm_genes_count <- length(dm_genes)

hyper_test <- function(input_table){
  dhyper(input_table["in_regulon", "dm"], dm_genes_count, 
         human_genes_count - dm_genes_count, sum(input_table["in_regulon", ]))
  
}

fisher_p <- function(input_table){
  fisher.test(input_table["in_regulon", "dm"], dm_genes_count, 
         human_genes_count - dm_genes_count, sum(input_table["in_regulon", ]))
  
}

fisher_p <- lapply(tables, fisher.test)

p_values <- lapply(tables, hyper_test)

adjusted_p <- p.adjust(p_values, "bonferroni")

hm_regulons <- names(adjusted_p[adjusted_p <= 0.1])

#----

gdata::keep(dm_genes, hm_regulons, adjusted_p, tables, tfs, sure = T)

save.image("./interanalysis_files/rdata_files/5_dm_regulons.RData")
