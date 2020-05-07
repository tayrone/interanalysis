library(gdata)

#subgroups <- c("wnt", "shh", "g3", "g4", "g34")

subgroups <- "g34"

load_mrs <- function(subgroup){
  
  load(paste0("./expression_analysis/rdata_files/network/", 
              subgroup, "_rtn.RData"))
  
  assign(paste0(subgroup, "_mrs"), 
         data.frame(tf = tna.get(rtna, what = "mra")[["Regulon"]], 
                    analysis = subgroup))
}

mrs_list <- lapply(subgroups, load_mrs)

mrs <- unlist(mrs_list)
mrs <- as.character(mrs)

gdata::keep(subgroups, mrs, sure = T)

#----

load("./methylation_analysis/control_g34_files/rdata_files/2_probewise_analysed.RData")
gdata::keep(genes_table, mrs, subgroups, sure = T)

dmgs <- names(genes_table)

any(duplicated(dmgs))

length(intersect(mrs, dmgs))

#----

load("/data4/tayrone25/expression_analysis/rdata_files/survival/g34_survival.RData")

gdata::keep(hazardous_regulons, dmgs, mrs, subgroups, sure = T)

sum(mrs %in% hazardous_regulons)
sum(hazardous_regulons %in% dmgs)

three_intersec <- intersect(intersect(mrs, hazardous_regulons), dmgs)

#---- Regulons that present coregulation ----

load("/data4/tayrone25/expression_analysis/rdata_files/duals/g34_duals.RData")

all(overlap$Regulon1 %in% overlap$Regulon2) #it's false, so we do the next line

coregulated_regulons <- unique(c(overlap$Regulon1, overlap$Regulon2))

intersect(coregulated_regulons, mrs)

intersect(coregulated_regulons, dmgs)

intersect(coregulated_regulons, hazardous_regulons)

intersect(three_intersec, coregulated_regulons)


hubs_coregulation <- sort(table(c(overlap$Regulon1, overlap$Regulon2)), decreasing = T)
sum(hubs_coregulation > 1)
