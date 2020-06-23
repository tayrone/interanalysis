library(RTN)
library(UpSetR)
library(dplyr)
library(vcd)
library(rcompanion)

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


#---- Then, check wich intersections are statistically significant,
# using the Cochran–Mantel–Haenszel test, which is basically chi-square
# for three dimensional tables ----

# This will make a contigency table that adresses all combinations 
# shown in the upsetr plot.
count_table <- data.frame(hm = c(1, 0, 1, 1, 0, 1, 0),
                          haz = c(0, 1, 1, 0, 0, 1, 1),
                          mr = c(rep(0, 3), rep(1, 4)))

colnames(count_table) <- c("hm", "haz", "mr")

complete_map_wide <- select(complete_map_wide, hm = `Highly Methylated Regulons`, 
                            haz = `Hazardous Regulons`, mr = `Master Regulators`)

count_table$count <- 
  sapply(rownames(count_table), 
        function(x) nrow(dplyr::inner_join(complete_map_wide, 
                                           count_table[x, ])))

table_mr <- xtabs(count ~ (hm + haz + mr), 
               data = count_table)
table_haz <- xtabs(count ~ (haz + mr + hm), 
                  data = count_table)
table_hm <- xtabs(count ~ (hm + mr + haz), 
                   data = count_table)

tables <- list(mr = table_mr, haz = table_haz, hm = table_hm)

results <- lapply(tables, mantelhaen.test)

#p-value is significant, so CMH tests are not really appropriate
woolf_results <- lapply(tables, woolf_test)

# for(i in 1:3){
#   print(groupwiseCMH(table_mr, group = i))
# }

gdata::keep(mrs, hazardous_regulons, hm_regulons, sure = T)

save.image("./interanalysis_files/rdata_files/7_hyper.RData")

