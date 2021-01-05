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

# another way to mapping the ensembl to the entrez
# kegg_dict <- bitr(AD$ensembl_gene_id, "ENSEMBL", "ENTREZID", org.Hs.eg.db)
# colnames(kegg_dict) <- c("ensembl_gene_id", "entrezgene_id")
# kegg_dict <- kegg_dict[!duplicated(kegg_dict$ensembl_gene_id),]
# kegg_dict <- kegg_dict[complete.cases(kegg_dict),]
# hgnc_symbol <- bitr(AD$ensembl_gene_id, "ENSEMBL", "SYMBOL", org.Hs.eg.db)
# colnames(hgnc_symbol) <- c("ensembl_gene_id", "hgnc_symbol")
# hgnc_symbol <- hgnc_symbol[!duplicated(hgnc_symbol$ensembl_gene_id),]
# hgnc_symbol <- hgnc_symbol[complete.cases(hgnc_symbol),]
# kegg_dict <- left_join(kegg_dict, hgnc_symbol, by = "ensembl_gene_id")

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
glyco_pathway <- enrich@result[grep("glycan", enrich@result$Description, ignore.case = TRUE),2:7] %>%
    tibble::rownames_to_column("kegg_id") 
    # mutate(Description = glue::glue("<a href='https://www.genome.jp/dbget-bin/www_bget?pathway:{kegg_id}'>{Description}</a>"))
glyco_set <- enrich@geneSets[rownames(enrich@result)[grep("glycan", enrich@result$Description, ignore.case = TRUE)]] 
pathway_id <- vector("integer", 0)
pathway_description <- vector("character", 0)
for (name in names(glyco_set)) {
    len = length(glyco_set[[name]])
    pathway = rep(name, len)
    pathway_id = c(pathway_id, pathway)
    description = glyco_pathway %>%
        dplyr::filter(kegg_id == name)
    description = description[["Description"]] %>%
        rep(len)
    pathway_description = c(pathway_description, description)
}
glyco_dict <- data.frame(
    entrezgene_id = as.integer(unlist(glyco_set)),
    pathway_id = pathway_id,
    pathway_description = pathway_description
) %>%
    dplyr::mutate(entrezgene_id = entrezgene_id) %>%
    left_join(kegg_dict, by = "entrezgene_id") %>%
    dplyr::filter(!is.na(ensembl_gene_id))
# ATTENTION: genes are involved in multiple pathways



# DE of all glycosylation related genes
AD_glyco <- AD %>%
    dplyr::filter(ensembl_gene_id %in% glyco_dict$ensembl_gene_id) %>%
    left_join(distinct(glyco_dict[,4:5]), by = "ensembl_gene_id")
# DE of glycosylation related genes with padj < 0.05
AD_glyco_filtered1 <- AD_glyco %>%
    filter(padj < 0.05) 
# DE of glycosylation related genes with padj < 0.05 and distinct to AD
distinctGenes <- DE[[region]]$distinctGenes
AD_glyco_filtered2 <- AD_glyco %>%
    filter(ensembl_gene_id %in% distinctGenes) %>%
    left_join(glyco_dict[,2:4], by = "ensembl_gene_id")
PA_glyco <- PA %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco$ensembl_gene_id)
PA_glyco_filtered1 <- PA %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
PA_glyco_filtered2 <- PA %>%
    dplyr::right_join(AD_glyco_filtered2[,1:2], by = "ensembl_gene_id")
PSP_glyco <- PSP %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco$ensembl_gene_id)
PSP_glyco_filtered1 <- PSP %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
PSP_glyco_filtered2 <- PSP %>%
    dplyr::right_join(AD_glyco_filtered2[,1:2], by = "ensembl_gene_id")

# CBE <- list(
#     AD_glyco = AD_glyco,
#     AD_glyco_sig = AD_glyco_filtered1,
#     AD_glyco_sig_exclu = AD_glyco_filtered2,
#     PA_glyco = PA_glyco,
#     PA_glyco_sig = PA_glyco_filtered1,
#     PA_glyco_sig_exclu = PA_glyco_filtered2,
#     PSP_glyco = PSP_glyco,
#     PSP_glyco_sig = PSP_glyco_filtered1,
#     PSP_glyco_sig_exclu = PSP_glyco_filtered2
# )
TCX <- list(
    AD_glyco = AD_glyco,
    AD_glyco_sig = AD_glyco_filtered1,
    AD_glyco_sig_exclu = AD_glyco_filtered2,
    PA_glyco = PA_glyco,
    PA_glyco_sig = PA_glyco_filtered1,
    PA_glyco_sig_exclu = PA_glyco_filtered2,
    PSP_glyco = PSP_glyco,
    PSP_glyco_sig = PSP_glyco_filtered1,
    PSP_glyco_sig_exclu = PSP_glyco_filtered2
)

# output <- list(CBE = CBE, TCX = TCX)
# output <- list(TCX = TCX)
# saveRDS(output, "../data/DE_glycosylation.rds")
# for (listname in names(output)) {
#     tmp <- output[[listname]]
#     for (tbl in names(tmp)) {
#         write.csv(tmp[[tbl]], paste0("../output/", paste(listname, tbl, sep = "_"), ".csv"))
#     }
# }
# DE_glycosylation <- readRDS("../data/DE_glycosylation.rds")

# 6 genes less
# old <- read.csv("../output/TCX_AD_glyco_sig_exclu.csv")
# dim(old)
# old$ensembl_gene_id[!(old$ensembl_gene_id %in% AD_glyco_filtered2$ensembl_gene_id)]
