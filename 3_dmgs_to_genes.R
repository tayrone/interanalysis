library(gdata)
library(dplyr)
library(biomaRt)
library(stringr)

load("/data4/tayrone25/methylation_analysis/control_g34_files/rdata_files/3_regionwise_analysed.RData")

#----

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

coord_genes <- 
  getBM(attributes = c("hgnc_symbol", "chromosome_name", 
                       "start_position", "end_position"),
      filters = "chromosomal_region", values = dmrs$coord, mart = grch37)

save(coord_genes, file = "3_biomart_result.RData")

#----

coord_genes$coord <- paste(coord_genes$chromosome_name, 
                           coord_genes$start_position,
                           coord_genes$end_position, sep = ":")

coord_genes <- dplyr::select(coord_genes, coord, hgnc_symbol)

test <- left_join(dmrs, coord_genes)







ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL")

filterlist <- list("6:33156164:33181870","protein_coding")

results <- getBM(attributes = c("hgnc_symbol","entrezgene", "chromosome_name",
                             "start_position", "end_position"),
              filters = "chromosomal_region", 
              values = "6:33156164:33181870", mart = ensembl54)






edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "chr6")
options(ucscChromosomeNames = T)



gnm <- GRanges("chr6:31618987-31639143")
gat <- GenomeAxisTrack(range = gnm)

gnm_gns <- getGeneRegionTrackForGviz(edbx, filter = GRangesFilter(gnm))

gtx <- GeneRegionTrack(gnm_gns, name = "tx", geneSymbol = TRUE,
                       showId = TRUE)

# load("/data4/tayrone25/methylation_analysis/control_g34_files/rdata_files/2_probewise_analysed.RData")
# 
# dm_cpgs <- dplyr::select(diff_methylated, islands_name, hgnc_symbol)
# 
# gdata::keep(dmrs, dm_cpgs, sure = T)


#----
# dmrs <- strsplit(dmrs$coord, "-", fixed = T)
# 
# dmrs <- data.frame(coord = sapply(dmrs, `[`, 1), stringsAsFactors = F)
# 
# dm_cpgs$islands_name <- strsplit(dm_cpgs$islands_name, "-", fixed = T)
# 
# dm_cpgs$islands_name <- sapply(dm_cpgs$islands_name, `[`, 1)
# 
# nrow(inner_join(dmrs, dm_cpgs, by = c("coord" = "islands_name")))


