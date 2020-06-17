library(RTN)
library(UpSetR)
library(dplyr)

load("./interanalysis_files/rdata_files/5_dm_regulons.RData")
load("../expression_analysis/rdata_files/survival/g34_survival.RData")
load("../expression_analysis/rdata_files/network/g34_rtn.RData")

#---- First step is to check intersections by regulons groups of interest ----

mrs <- tna.get(rtna, what = "mra")

gdata::keep(mrs, hazardous_regulons, hm_regulons, sure = T)

long_data <- 
  data.frame(tf = c(mrs$Regulon, hazardous_regulons, hm_regulons),
             analysis = c(rep("Master Regulators", length(mrs$Regulon)),
                          rep("Hazardous Regulons", length(hazardous_regulons)),
                          rep("Highly Methylated Regulons", length(hm_regulons))))

complete_map_wide <- dcast(long_data, tf~analysis, length, 
                           value.var = "analysis")   

rownames(complete_map_wide) <- complete_map_wide$tf
complete_map_wide <- select(complete_map_wide, -tf)

upset(complete_map_wide, order.by = "freq")


#---- Then, check wich intersections are statistically significant ----

# Tip: check for 3 variable chi square test

count_table <- data.frame(matrix(0, nrow = 2, ncol = 2))

count_table[1, 1] <- sum((complete_map_wide$`Hazardous Regulons` + 
                          complete_map_wide$`Highly Methylated Regulons`) == 2 &
                          complete_map_wide$`Master Regulators` == 0) 
                     
colnames(count_table) <- c("is_hm", "not_hm")
rownames(count_tables) <- c("is_haz", "not_haz")


