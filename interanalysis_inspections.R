library(gdata)

subgroup <- "g34" #analysis under inspection


#---- Expression analysis results ----

load(paste0("./expression_analysis/rdata_files/network/", 
            subgroup, "_rtn.RData"))
  
mrs <-  tna.get(rtna, what = "mra")[["Regulon"]]

gdata::keep(mrs, sure = T)


#---- Probewise methylation analysis results ----

load("./methylation_analysis/control_g34_files/rdata_files/2_probewise_analysed.RData")
gdata::keep(genes_table, mrs, sure = T)

dmgs <- names(genes_table)

any(duplicated(dmgs))

length(intersect(mrs, dmgs))


#---- Survival analysis results ----

load("./expression_analysis/rdata_files/survival/g34_survival.RData")

gdata::keep(hazardous_regulons, dmgs, mrs, sure = T)

sum(mrs %in% hazardous_regulons)
sum(hazardous_regulons %in% dmgs)

three_intersec <- intersect(intersect(mrs, hazardous_regulons), dmgs)


#---- Coregulation analysis results ----

load("./expression_analysis/rdata_files/duals/g34_duals.RData")

gdata::keep(overlap, dmgs, hazardous_regulons, mrs, three_intersec, sure = T)

all(overlap$Regulon1 %in% overlap$Regulon2) #it's false, so we do the next line

coregulated_regulons <- unique(c(overlap$Regulon1, overlap$Regulon2))

intersect(coregulated_regulons, mrs)

intersect(coregulated_regulons, dmgs)

intersect(coregulated_regulons, hazardous_regulons)

intersect(three_intersec, coregulated_regulons)


hubs_coregulation <- sort(table(c(overlap$Regulon1, overlap$Regulon2)), decreasing = T)
sum(hubs_coregulation > 1)
