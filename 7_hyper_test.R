library(hypeR)

load("../expression_analysis/rdata_files/tfs.RData")

symbols <- list(all = tfs)

regulons_list <- list(mrs = mrs$Regulon, haz = hazardous_regulons, 
                      hm = hm_regulons)

hyp <- hypeR(symbols, regulons_list, 1639, verbose = T)
