#----

dat <- data.frame(
  rm = sum(complete_map_wide$reguladores_mestres == 1 & 
             complete_map_wide$regulons_de_risco == 1),
  not_rm = sum(complete_map_wide$reguladores_mestres == 0 & 
                 complete_map_wide$regulons_de_risco == 1)
)
dat <- rbind(dat, 
             c(sum(complete_map_wide$reguladores_mestres == 1 & 
                     complete_map_wide$regulons_de_risco == 0),
               sum(complete_map_wide$reguladores_mestres == 0 & 
                     complete_map_wide$regulons_de_risco == 0)))
rownames(dat) <- c("de_risco", "not_de_risco")

chisq.test(dat)$expected
