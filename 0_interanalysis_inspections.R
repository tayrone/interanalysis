library(gdata)

#---- To make this script analysis-independent ----

expression_subgroup <- "g34" 
methylation_subgroup <- "control_g34"

#---- Expression analysis results ----

load(paste0("./expression_analysis/rdata_files/network/", 
            expression_subgroup, "_rtn.RData"))
  
mrs <-  tna.get(rtna, what = "mra")[["Regulon"]]

gdata::keep(mrs, expression_subgroup, methylation_subgroup, sure = T)


#---- Probewise methylation analysis results ----

load(paste0("./methylation_analysis/", methylation_subgroup,
            "_files/rdata_files/2_probewise_analysed.RData"))

gdata::keep(genes_table, mrs, expression_subgroup, methylation_subgroup, sure = T)

dmgs <- names(genes_table)

any(duplicated(dmgs))

length(intersect(mrs, dmgs))


#---- Survival analysis results ----

load(paste0("./expression_analysis/rdata_files/survival/", expression_subgroup,
            "_survival.RData"))

gdata::keep(hazardous_regulons, dmgs, mrs, expression_subgroup, 
            methylation_subgroup, sure = T)

sum(mrs %in% hazardous_regulons)
sum(hazardous_regulons %in% dmgs)

three_intersec <- intersect(intersect(mrs, hazardous_regulons), dmgs)


#---- Coregulation analysis results ----

load(paste0("./expression_analysis/rdata_files/duals/", expression_subgroup,
            "_duals.RData"))

gdata::keep(overlap, dmgs, hazardous_regulons, mrs, three_intersec, sure = T)

all(overlap$Regulon1 %in% overlap$Regulon2) #it's false, so we do the next line

coregulated_regulons <- unique(c(overlap$Regulon1, overlap$Regulon2))

intersect(coregulated_regulons, mrs)

intersect(coregulated_regulons, dmgs)

intersect(coregulated_regulons, hazardous_regulons)

intersect(three_intersec, coregulated_regulons)


hubs_coregulation <- sort(table(c(overlap$Regulon1, overlap$Regulon2)), 
                          decreasing = T)
sum(hubs_coregulation > 1)


save.image("./interanalysis_files/rdata_files/0_inspections_image.RData")
