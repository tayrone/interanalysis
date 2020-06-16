# This script aims to identify highly methylated regulons among
# the ones identified as master regulators.

options(stringsAsFactors = F)
library(tidyverse)
library(gdata)
library(ggplot2)

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

sig_test <- function(x){
  
  x <- as.data.frame(x)
  
  if(any(chisq.test(x)$expected < 5)){
    return(fisher.test(x)$p.value)
  }else{
    return(chisq.test(x)$p.value)
  }
  
}

p_values <- sapply(tables, sig_test)

adjusted_p <- p.adjust(p_values, "bonferroni")

length(unique(c(tfs, dm_genes)))

dm_regulons <- names(adjusted_p[adjusted_p <= 0.1])

#---- Check if more than 50% of regulons' elements are DM ----

regulon_methylation <- function(x, threshold){ 
  x <- as.data.frame(x)
  regulon_elements <- x["in_regulon", ]
  
  if((regulon_elements$dm/(regulon_elements$not_dm + regulon_elements$dm)) > threshold &
     (regulon_elements$dm/(regulon_elements$not_dm + regulon_elements$dm)) <= (threshold + 0.05)){

    return (TRUE)
  }else{
    return (FALSE)
  }
}

plot_data <- data.frame(threshold = NULL, hm_regulons = NULL)


threshold <- seq(0, 0.95, 0.05)

for(i in threshold){
  
  hm_regulons <- sapply(tables, regulon_methylation, threshold = i)
  hm_regulons <- names(tables)[hm_regulons]
  
  current_values <- 
    data.frame(threshold = i, 
               hm_regulons = length(intersect(hm_regulons, dm_regulons)))
  
  plot_data <- rbind(plot_data, current_values)
}

ggplot(plot_data, aes(x = threshold, y = hm_regulons)) +
  geom_col(orientation = "x", position = position_nudge(x = 0.025)) +
  labs(x = "Proporção de elementos diferencialmente metilados, no regulon",
       y = "Número de regulons") +
  theme_minimal()


#----

gdata::keep(dmrs, gene_coords, dm_regulons, p_values, sure = T)

save.image("./interanalysis_files/rdata_files/5_dm_regulons.RData")
