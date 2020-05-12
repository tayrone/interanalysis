library(gdata)

load("./plots_image.RData")

gdata::keep(complete_map_wide, sure = T)

colnames(complete_map_wide) <- c("mrs", "dmgs", "haz", "coreg")


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

rm(complete_map_wide)


#---- The next loop is necessary to run the significance test for all
# generated data frames ----

for(object in ls()){
  
  cat(as.name(object), "\nExpected: \n")
  
  print(chisq.test(get(object))$expected)
  
  cat("\nChi test: \n")
  
  print(chisq.test(get(object)))
  
  cat("\nP-value: \n")
  
  print(chisq.test(get(object))$p.value)
  
  cat("\n----------\n")
}

rm(object)
