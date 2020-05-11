#---- This script build a set intersection graph. It takes as
# input all transcription factors identified as master regulators,
# for each one of the analyses. -----

library(reshape2)
library(UpSetR)
library(tidyverse)


load("inspections_image.RData")

reg_elements <- data.frame(element = mrs, 
                           analysis = "reguladores_mestres")
reg_elements <- rbind(reg_elements, data.frame(element = dmgs, 
                                               analysis = "diferencialmente_metilados"))
reg_elements <- rbind(reg_elements, data.frame(element = hazardous_regulons, 
                                               analysis = "regulons_de_risco"))
reg_elements <- rbind(reg_elements, data.frame(element = coregulated_regulons, 
                                               analysis = "correguladores"))


#---- Wide format is required by upset method ----

complete_map_wide <- dcast(reg_elements, element~analysis, length, 
                           value.var = "analysis")

rownames(complete_map_wide) <- complete_map_wide$element
complete_map_wide <- select(complete_map_wide, -element)


upset(complete_map_wide, order.by = "freq")


