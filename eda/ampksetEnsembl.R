setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("dplyr","HTSet", "clusterProfiler", "biomaRt", "org.Hs.eg.db")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

DE <- readRDS("../data/DE.rds")
# run once----
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#-------------

region <- "TCX"
AD <- DE[[region]]$AD
PA <- DE[[region]]$PA
PSP <- DE[[region]]$PSP

kegg_dict <- getBM(
    filters="ensembl_gene_id",
    attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
    values=AD$ensembl_gene_id,
    mart=mart)
kegg_dict <- kegg_dict[!duplicated(kegg_dict$ensembl_gene_id),]
kegg_dict <- kegg_dict[complete.cases(kegg_dict),]

getGeneList <- function(DETable = AD){
    # return a geneList with decreasing logFC
    DETable <- dplyr::filter(DETable, padj<0.05)
    table2 <- data.frame(
        ensembl_gene_id = DETable$ensembl_gene_id,
        logFC = DETable$logFC
    ) %>%
        right_join(kegg_dict, by = "ensembl_gene_id")
    geneList <- table2$logFC
    names(geneList) <- table2$entrezgene_id
    geneList <- sort(geneList, decreasing = T)
}
geneList <- getGeneList(AD)




enrich <- enrichKEGG(
    gene = names(geneList),
    organism = 'hsa'
)
ampk_set <- enrich@geneSets$hsa04152



ampk_dict <- data.frame(
    entrezgene_id = as.integer(ampk_set),
    pathway_id = "hsa04152",
    pathway_description = enrich@result["hsa04152","Description"]
) %>%
    # dplyr::mutate(entrezgene_id = entrezgene_id) %>%
    left_join(kegg_dict, by = "entrezgene_id") %>%
    dplyr::filter(!is.na(ensembl_gene_id))


AD_ampk <- AD %>%
    dplyr::filter(ensembl_gene_id %in% ampk_dict$ensembl_gene_id) %>%
    left_join(distinct(ampk_dict[,4:5]), by = "ensembl_gene_id")
# DE of glycosylation related genes with padj < 0.05
AD_ampk_filtered1 <- AD_ampk %>%
    filter(padj < 0.05) 
# DE of glycosylation related genes with padj < 0.05 and distinct to AD
distinctGenes <- DE[[region]]$distinctGenes
AD_ampk_filtered2 <- AD_ampk %>%
    filter(ensembl_gene_id %in% distinctGenes) %>%
    left_join(ampk_dict[,2:4], by = "ensembl_gene_id")
PA_ampk <- PA %>%
    dplyr::filter(ensembl_gene_id %in% AD_ampk$ensembl_gene_id)
PA_ampk_filtered1 <- PA %>%
    dplyr::filter(ensembl_gene_id %in% AD_ampk_filtered1$ensembl_gene_id)
PA_ampk_filtered2 <- PA %>%
    dplyr::right_join(AD_ampk_filtered2[,1:2], by = "ensembl_gene_id")
PSP_ampk <- PSP %>%
    dplyr::filter(ensembl_gene_id %in% AD_ampk$ensembl_gene_id)
PSP_ampk_filtered1 <- PSP %>%
    dplyr::filter(ensembl_gene_id %in% AD_ampk_filtered1$ensembl_gene_id)
PSP_ampk_filtered2 <- PSP %>%
    dplyr::right_join(AD_ampk_filtered2[,1:2], by = "ensembl_gene_id")

TCX <- list(
    AD_ampk = AD_ampk,
    AD_ampk_sig = AD_ampk_filtered1,
    AD_ampk_sig_exclu = AD_ampk_filtered2,
    PA_ampk = PA_ampk,
    PA_ampk_sig = PA_ampk_filtered1,
    PA_ampk_sig_exclu = PA_ampk_filtered2,
    PSP_ampk = PSP_ampk,
    PSP_ampk_sig = PSP_ampk_filtered1,
    PSP_ampk_sig_exclu = PSP_ampk_filtered2
)

# output <- list(CBE = CBE, TCX = TCX)
output <- list(TCX = TCX)
saveRDS(output, "../data/DE_ampk.rds")
for (listname in names(output)) {
    tmp <- output[[listname]]
    for (tbl in names(tmp)) {
        write.csv(tmp[[tbl]], paste0("../output/", paste(listname, tbl, sep = "_"), ".csv"))
    }
}



