#---- This script build a set intersection graph. It takes as
# input all transcription factors identified as master regulators,
# for each one of the analyses. -----

library(reshape2)
library(UpSetR)
library(tidyverse)
library(RTN)

#subgroups <- c("wnt", "shh", "g3", "g4", "g34")
subgroups <- "g34"

load_mrs <- function(subgroup){
  
  load(paste0("../expression_analysis/rdata_files/network/", 
              subgroup, "_rtn.RData"))
  
  assign(paste0(subgroup, "_mrs"), 
         data.frame(tf = tna.get(rtna, what = "mra")[["Regulon"]], 
                    analysis = subgroup))
}

mrs_list <- lapply(subgroups, load_mrs)

rm(rtna, rtni, subgroups, load_mrs)

#---- Organize all data on a list ----

load(paste0("../expression_analysis/rdata_files/survival/", 
            subgroup, "_rtn.RData"))

load(paste0("../methylation_analysis/control_", 
            subgroup, "_files/rdata_files/", "2_probewise_analysed.RData"))

gdata::keep(diff_methylated, mrs_list, hazardous_regulons, subgroups, sure = T)


dm_genes <- unique(diff_methylated$hgnc_symbol)

mrs_list[[1]]$analysis <- "Master Regulators"

mrs_list[[2]] <- data.frame(tf = hazardous_regulons, analysis = "Hazardous Regulons")
mrs_list[[3]] <- data.frame(tf = dm_genes, analysis = "Differentially Methylated Genes")

complete_map <- do.call("rbind", mrs_list)


#---- Wide format is required by upset method ----

complete_map_wide <- dcast(complete_map, tf~analysis, length, 
                           value.var = "analysis")

# only one observation is dropped out :)
complete_map_wide <- complete_map_wide[complete.cases(complete_map_wide), ]

rownames(complete_map_wide) <- complete_map_wide$tf
complete_map_wide <- select(complete_map_wide, -tf)


upset(complete_map_wide, order.by = "freq")
 
