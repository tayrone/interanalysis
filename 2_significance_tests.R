library(gdata)

#---- To make this script analysis-independent ----

expression_subgroup <- "g34" 

load("./interanalysis_files/rdata_files/1_plots_image.RData")
load(paste0("../expression_analysis/rdata_files/network/", expression_subgroup,
            "_rtn.RData"))

gdata::keep(complete_map_wide, rtni, sure = T)

colnames(complete_map_wide) <- c("mrs", "dmgs", "haz", "coreg")


#---- TFs which were left unnoticed by all four analyses are binded
# to the data frame, containing all 0 rows, since they might 
# influence on tests results----

tfs <- rtni@regulatoryElements
names(tfs) <- NULL

unnoticed_tfs <- tfs[!(tfs %in% rownames(complete_map_wide))]

unnoticed_tfs <- data.frame(matrix(0, nrow = length(unnoticed_tfs), 
                                   ncol = length(colnames(complete_map_wide))),
                            row.names = unnoticed_tfs)

colnames(unnoticed_tfs) <- colnames(complete_map_wide)

complete_map_wide <- rbind(complete_map_wide, unnoticed_tfs)

rm(rtni, tfs)


#---- The following nested loops are required to build contingency tables
# for each possible pair for all complete_map_wide variables ----

for(i in 1:ncol(complete_map_wide)){
  
  for(j in 1:ncol(complete_map_wide)){
    
    if((i != j) & (i < j)){
      
      freqs <- data.frame(matrix(NA_integer_, nrow = 2, ncol = 2))
      
      colnames(freqs) <- c(paste0(colnames(complete_map_wide)[i]),
                           paste0("not_", colnames(complete_map_wide)[i]))
      
      rownames(freqs) <- c(paste0(colnames(complete_map_wide)[j]),
                           paste0("not_", colnames(complete_map_wide)[j]))
      
      
      freqs[1, 1] <- sum(complete_map_wide[, i] == 1 &
                           complete_map_wide[, j] == 1)
      
      freqs[2, 1] <- sum(complete_map_wide[, i] == 1 &
                           complete_map_wide[, j] == 0)
      
      freqs[1, 2] <- sum(complete_map_wide[, i] == 0 &
                           complete_map_wide[, j] == 1)
      
      freqs[2, 2] <- sum(complete_map_wide[, i] == 0 &
                           complete_map_wide[, j] == 0)
      
      assign(paste(colnames(complete_map_wide)[i], 
                   colnames(complete_map_wide)[j], sep = "_"),
             freqs)
      
      rm(freqs)
      
    }
  }
}

rm(complete_map_wide, i, j, unnoticed_tfs)


#---- The next loop is necessary to run a significance test for all
# generated data frames. If any cell of the Chi-square expected 
# values data frame is less than 5, then Fisher's exact test is 
# more appropriate ----
# 
# for(object in ls()){
# 
#   cat(as.name(object), "Expected: \n")
# 
#   print(chisq.test(get(object))$expected)
# 
#   if(any(chisq.test(get(object))$expected < 5)){
# 
#     cat("Fisher test: \n")
#     print(fisher.test(get(object)))
# 
#     cat("P-value: \n")
#     print(fisher.test(get(object))$p.value)
# 
#   }else{
# 
#     cat("Chi test: \n")
#     print(chisq.test(get(object)))
# 
#     cat("P-value: \n")
#     print(chisq.test(get(object))$p.value)
# 
#   }
#   cat("\n----------\n")
# }


rm(object)
