library(ggplot2)
load("./interanalysis_files/rdata_files/5_dm_regulons.RData")

regulon_methylation <- function(x, threshold){ 
  x <- as.data.frame(x)
  regulon_elements <- x["in_regulon", ]
  
  if((regulon_elements$dm/(regulon_elements$not_dm + regulon_elements$dm)) > threshold &
     (regulon_elements$dm/(regulon_elements$not_dm + regulon_elements$dm)) <= (threshold + 0.05)){
    
    return (TRUE)
  }else{
    return (FALSE)
  }
}

plot_data <- data.frame(threshold = NULL, hm_regulons = NULL)


threshold <- seq(0, 0.95, 0.05)

for(i in threshold){
  
  hm_regulons <- sapply(tables, regulon_methylation, threshold = i)
  hm_regulons <- names(tables)[hm_regulons]
  
  current_values <- 
    data.frame(threshold = i, 
               hm_regulons = length(intersect(hm_regulons, dm_regulons)))
  
  plot_data <- rbind(plot_data, current_values)
}

ggplot(plot_data, aes(x = threshold, y = hm_regulons)) +
  geom_col(orientation = "x", position = position_nudge(x = 0.025)) +
  labs(x = "Proporção de elementos diferencialmente metilados, no regulon",
       y = "Número de regulons") +
  theme_minimal()
