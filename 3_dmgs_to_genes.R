# The following code maps differentially methylated regions to its genes.

library(gdata)
library(dplyr)
library(biomaRt)
library(stringr)

load("./control_g34_files/rdata_files/3_regionwise_analysed.RData")

#---- A specific way to write chromosome coordinates is required by BioMart ----

dmrs <- dplyr::select(dm_regions$results, coord, no.cpgs)
gdata::keep(dmrs, sure = T)

dmrs$coord <- str_replace(dmrs$coord, "chr", "")
dmrs$coord <- str_replace(dmrs$coord, "-", ":")

rownames(dmrs) <- NULL

#---- Makes sure hg19/grch37 genome is being used, like all other
# coordinate objects from methylation analysis ----

grch37 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                 host = "grch37.ensembl.org", path = "/biomart/martservice", 
                 dataset = "hsapiens_gene_ensembl")

dmr_symbol <- NULL
gene_coords <- NULL


#---- This for loop is applied so we know which gene symbols belong
# to each DMR, exactly ----

for(x in dmrs$coord){

  bm_out <- 
    getBM(attributes = c("hgnc_symbol", "chromosome_name", 
                         "start_position", "end_position"),
          filters = "chromosomal_region", values = x, mart = grch37)
  
  gene_coords <- rbind(gene_coords, bm_out)
  
  dmr_symbol <- append(dmr_symbol, str_c(bm_out$hgnc_symbol, collapse = " "))

}


get_genes <- function(x){
  bm_out <- 
    getBM(attributes = c("hgnc_symbol", "chromosome_name", 
                         "start_position", "end_position"),
          filters = "chromosomal_region", values = x, mart = grch37)
}


gene_coords <- lapply(dmrs$coord, get_genes)

gene_coords <- lapply(gene_coords, function(x) mutate_all(x, as.character))

gene_coords <- dplyr::bind_rows(gene_coords, .id = "dmr_id")

dmrs$dmr_id <- rownames(dmrs)



save(gene_coords, dmrs, 
     file = "./interanalysis_files/rdata_files/3_dmgs_to_genes.RData")

