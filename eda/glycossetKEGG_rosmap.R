setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("dplyr", "HTSet", "clusterProfiler", "biomaRt", "org.Hs.eg.db", "xlsx")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

DE <- readRDS("../data/DE_rosmap.rds")
load("glycanSet.rda")
# run once----
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#-------------

region <- "dorsolateral"
AD <- DE[[region]]$AD
MCI <- DE[[region]]$MCI
Other <- DE[[region]]$Other
anova <- DE[[region]]$anova

# The DE data has ensembl id and hgnc symbol. Now find entrez gene id for these
# entries: 
# Why need entrez gene id? KEGG da is based on entrez gene id.
ensembl_to_entrez <- getBM(
    filters="ensembl_gene_id",
    attributes=c("ensembl_gene_id", "entrezgene_id"),
    values=AD$ensembl_gene_id,
    mart=mart)
ensembl_to_entrez <- ensembl_to_entrez[!duplicated(ensembl_to_entrez$ensembl_gene_id),]
ensembl_to_entrez <- ensembl_to_entrez[complete.cases(ensembl_to_entrez),]

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

#### Now we got entrez id, pathway id, and pathway description. 
#### Map this information to DE table. First, we need to convert entrez id 
#### to ensembl id:
entrez_to_ensembl <- data.frame(
    entrezgene_id = as.integer(unlist(glyco_set)),
    pathway_id = pathway_id,
    pathway_description = pathway_description
) %>%
    left_join(ensembl_to_entrez, by = "entrezgene_id") %>%
    dplyr::filter(!is.na(ensembl_gene_id))
# ATTENTION: genes are involved in multiple pathways

# DE of all glycosylation related genes
AD_glyco <- AD %>%
    dplyr::filter(ensembl_gene_id %in% entrez_to_ensembl$ensembl_gene_id) %>%
    right_join(distinct(entrez_to_ensembl[,2:4]), by = "ensembl_gene_id") %>%
    dplyr::select(1:4, 9:10, 5:8)
# DE of glycosylation related genes with padj < 0.05
AD_glyco_filtered1 <- AD_glyco %>%
    dplyr::filter(padj < 0.05) 
# DE of glycosylation related genes with padj < 0.05 and distinct to AD
distinctGenes <- DE[[region]]$distinctGenes
AD_glyco_filtered2 <- AD_glyco %>%
    filter(ensembl_gene_id %in% distinctGenes)

MCI_glyco <- MCI %>%
    right_join(AD_glyco[,c(1,2,5,6)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
    dplyr::select(1:4, 9:10, 5:8)
MCI_glyco_filtered1 <- MCI_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
MCI_glyco_filtered2 <- MCI_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered2$ensembl_gene_id)

Other_glyco <- Other %>%
    right_join(AD_glyco[,c(1,2,5,6)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
    dplyr::select(1:4, 9:10, 5:8)
Other_glyco_filtered1 <- Other_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
Other_glyco_filtered2 <- Other_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered2$ensembl_gene_id)

anova_glyco <- anova %>%
    right_join(AD_glyco[,c(1,2,5,6)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
    dplyr::select(1:4, 8:9, 5:7)
anova_glyco_filtered1 <- anova_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
anova_glyco_filtered2 <- anova_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered2$ensembl_gene_id)


bm <- list(
    AD_glyco = AD_glyco,
    AD_glyco_sig = AD_glyco_filtered1,
    AD_glyco_sig_exclu = AD_glyco_filtered2,
    MCI_glyco = MCI_glyco,
    MCI_glyco_sig = MCI_glyco_filtered1,
    MCI_glyco_sig_exclu = MCI_glyco_filtered2,
    Other_glyco = Other_glyco,
    Other_glyco_sig = Other_glyco_filtered1,
    Other_glyco_sig_exclu = Other_glyco_filtered2,
    anova_glyco = anova_glyco,
    anova_glyco_sig = anova_glyco_filtered1,
    anova_glyco_sig_exclu = anova_glyco_filtered2
)

wb <- createWorkbook()
sheet <- createSheet(wb, "glyco_sig_exclu")
addDataFrame(bm$AD_glyco_sig_exclu, sheet=sheet, startColumn=1, row.names=T)
addDataFrame(bm$MCI_glyco_sig_exclu[,7:10], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2, row.names=FALSE)
addDataFrame(bm$Other_glyco_sig_exclu[,7:10], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2+4, row.names=FALSE)
addDataFrame(bm$anova_glyco_sig_exclu[,7:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2+4+4, row.names=FALSE)

sheet <- createSheet(wb, "glyco_sig")
addDataFrame(bm$AD_glyco_sig, sheet=sheet, startColumn=1, row.names=T)
addDataFrame(bm$MCI_glyco_sig[, 7:10], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig)+2, row.names=FALSE)
addDataFrame(bm$Other_glyco_sig[, 7:10], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig)+2+4, row.names=FALSE)
addDataFrame(bm$anova_glyco_sig[,7:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig)+2+4+4, row.names=FALSE)

sheet <- createSheet(wb, "glyco")
addDataFrame(bm$AD_glyco, sheet=sheet, startColumn=1, row.names=T)
addDataFrame(bm$MCI_glyco[,7:10], sheet=sheet, startColumn=ncol(bm$AD_glyco)+2, row.names=FALSE)
addDataFrame(bm$Other_glyco[,7:10], sheet=sheet, startColumn=ncol(bm$AD_glyco)+2+4, row.names=FALSE)
addDataFrame(bm$anova_glyco[,7:9], sheet=sheet, startColumn=ncol(bm$AD_glyco)+2+4+4, row.names=FALSE)

saveWorkbook(wb, paste0("../output/rosmap_modified_", region, ".xlsx"))