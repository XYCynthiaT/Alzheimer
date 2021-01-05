library(clusterProfiler)
library(biomaRt)
library(dplyr)

genes <- readRDS("../data/DE.rds")[["TCX"]]$AD

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_to_entrez <- getBM(
    filters="ensembl_gene_id",
    attributes=c("ensembl_gene_id", "entrezgene_id"),
    values=genes$ensembl_gene_id,
    mart=mart)
ensembl_to_entrez <- ensembl_to_entrez[!duplicated(ensembl_to_entrez$ensembl_gene_id),]
ensembl_to_entrez <- ensembl_to_entrez[complete.cases(ensembl_to_entrez),]

## Pathway Enrichment Analysis------------------------
getGeneList <- function(DETable = genes){
    # return a geneList with decreasing logFC
    # DETable <- dplyr::filter(DETable, padj<0.05)
    table2 <- data.frame(
        ensembl_gene_id = DETable$ensembl_gene_id,
        logFC = DETable$logFC
    ) %>%
        right_join(ensembl_to_entrez, by = "ensembl_gene_id")
    geneList <- table2$logFC
    names(geneList) <- table2$entrezgene_id
    geneList <- sort(geneList, decreasing = T)
}
geneList <- getGeneList()

enrich <- enrichKEGG(
    gene = names(geneList),
    organism = 'hsa'
)

## Extract entrez gene ids of glycan related genes from enrichment results:
## Improve it##
glyco_description <- grep("glycan", enrich@result$Description, ignore.case = TRUE)
glyco_set <- enrich@geneSets[rownames(enrich@result)[glyco_description]]

pathway_id <- vector("integer", 0)
pathway_description <- vector("character", 0)
for (name in names(glyco_set)) {
    len = length(glyco_set[[name]])
    pathway = rep(name, len)
    pathway_id = c(pathway_id, pathway)
    description = enrich@result$Description[rownames(enrich@result) == name] %>%
        rep(len)
    pathway_description = c(pathway_description, description)
}

save(glyco_set, pathway_id, pathway_description, file = "glycanSet.rda")
