library(RTN)
library(UpSetR)
library(dplyr)
library(SuperExactTest)
library(reshape2)

load("./interanalysis_files/rdata_files/5_dm_regulons.RData")
load("../expression_analysis/rdata_files/survival/g34_survival.RData")
load("../expression_analysis/rdata_files/network/g34_rtn.RData")

#---- First step is to check intersections by regulons groups of interest ----

mrs <- tna.get(rtna, what = "mra")

gdata::keep(mrs, hazardous_regulons, hm_regulons, tfs, sure = T)

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

#---- Peforms the hypergeometric test for defining statistical significance of 
# set overlapping groups ----

input <- list(mrs = mrs$Regulon, hm = hm_regulons, haz = hazardous_regulons)

result <- supertest(input, n = length(tfs))

adjusted_intersections_p <- p.adjust(result$P.value, "bonferroni")

# plot(result, Layout = "landscape", sort.by = "size", keep=FALSE,
#      bar.split = c(70,180), show.elements = F, elements.cex = 0.7,
#      elements.list = subset(summary(result)$Table, Observed.Overlap <= 10),
#      show.expected.overlap = TRUE, expected.overlap.style = "hatchedBox",
#      color.expected.overlap = 'red')

gdata::keep(adjusted_intersections_p, mrs, result, hazardous_regulons, 
            hm_regulons, sure = T)

save.image("./interanalysis_files/rdata_files/6_mr_intersections.RData")
