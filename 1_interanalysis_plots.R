# This script build a set intersection graph. It takes as
# input all transcription factors identified as master regulators,
# for each one of the analyses. 

library(reshape2)
library(UpSetR)
library(tidyverse)

load("./interanalysis_files/rdata_files/0_inspections_image.RData")

reg_elements <- rbind(data.frame(element = mrs, 
                                 analysis = "Master Regulators"), 
                      data.frame(element = dmgs, 
                                 analysis = "Diferentially Methylated Elements"),
                      data.frame(element = hazardous_regulons, 
                                 analysis = "Hazardous Regulons"),
                      data.frame(element = coregulated_regulons, 
                                 analysis = "Co-regulating Elements"))


#---- Wide format is required by upset method ----

complete_map_wide <- dcast(reg_elements, element~analysis, length, 
                           value.var = "analysis")

rownames(complete_map_wide) <- complete_map_wide$element
complete_map_wide <- select(complete_map_wide, -element)

png("./interanalysis_files/plots/upset_plot.png", units = "in", width = 9, height = 5, res = 300)

upset(complete_map_wide, order.by = "freq")

graphics.off()


save.image("./interanalysis_files/rdata_files/1_plots_image.RData")

