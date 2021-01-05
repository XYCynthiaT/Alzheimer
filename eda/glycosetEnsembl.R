setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("dplyr","HTSet", "clusterProfiler", "biomaRt")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

ensembl <- read.csv("../raw-data/ensemblOutput.csv") %>%
    dplyr::rename(ensembl_gene_id = id_with_url)
DE <- readRDS("../data/DE.rds")
ADfiltered <- DE$CBE$ADflitered
# run once----
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#-------------

ids <- getBM(
    filters="ensembl_gene_id",
    attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
    values=ensembl$ensembl_gene_id,
    mart=mart)
ids <- ids[!duplicated(ids$ensembl_gene_id),]
table <- merge(ensembl, ids, by = "ensembl_gene_id")
# ensembl table contains several novel proteins/genes 
# that don't have corresponding entrez_id

AD_glyco_filterd <- ADfiltered %>%
    filter(ensembl_gene_id %in% table$ensembl_gene_id) %>%
    left_join(table, by = "ensembl_gene_id") %>%
    dplyr::select(1:3, 8:10, 4:7)

# getGeneList <- function(DETable = ADfiltered){
#     # return a geneList with decreasing logFC
#     table2 <- data.frame(
#         ensembl_gene_id = DETable$ensembl_gene_id,
#         logFC = ADfiltered$logFC
#     ) %>%
#         right_join(table, by = "ensembl_gene_id")
#     geneList <- table2$logFC
#     names(geneList) <- table2$entrezgene_id
#     geneList <- sort(geneList, decreasing = T)
# }
# geneList <- getGeneList(ADfiltered)
#  
# glycosylation <- gseKEGG(
#     geneList = geneList,
#     organism = 'hsa',
#     pvalueCutoff = 1
# )
# tmp <- glycosylation@result[grep("glycan", glycosylation@result$Description, ignore.case = TRUE),2:7] %>%
#     tibble::rownames_to_column("kegg_id") %>%
#     mutate(Description = glue::glue("<a href='https://www.genome.jp/dbget-bin/www_bget?pathway:{kegg_id}'>{Description}</a>"))

# library("pathview")
# pathview(gene.data  = geneList,
#                      pathway.id = "hsa00514",
#                      species    = "hsa")