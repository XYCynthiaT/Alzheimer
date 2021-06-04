setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("dplyr", "HTSet", "xlsx")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

DE <- readRDS("../data/DE_rosmap_adj.rds")
# reactome <- readRDS("../data/glyco_reactome.rds") %>%
#     dplyr::rename(ensembl_gene_id = Ensembl)
reactome <- readRDS("../data/glyco.rds") %>%
    dplyr::rename(ensembl_gene_id = Ensembl)


region <- "dorsolateral"
AD <- DE[[region]]$AD
MCI <- DE[[region]]$MCI
Other <- DE[[region]]$Other
anova <- DE[[region]]$anova

# DE of all glycosylation related genes
AD_glyco <- AD %>%
    dplyr::filter(ensembl_gene_id %in% reactome$ensembl_gene_id) %>%
    left_join(distinct(reactome[,c(1:3, 6)]), by = "ensembl_gene_id") %>%
    dplyr::select(1:2, 4, 9:11, 5:8)
# DE of glycosylation related genes with padj < 0.05
AD_glyco_filtered1 <- AD_glyco %>%
    dplyr::filter(padj < 0.05) 
# DE of glycosylation related genes with padj < 0.05 and distinct to AD
distinctGenes <- DE[[region]]$distinctGenes
AD_glyco_filtered2 <- AD_glyco %>%
    dplyr::filter(ensembl_gene_id %in% distinctGenes)

MCI_glyco <- MCI %>%
    right_join(AD_glyco[,c(1,2,4:5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
    dplyr::select(1:2, 4, 9:10, 5:8)
MCI_glyco_filtered1 <- MCI_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
MCI_glyco_filtered2 <- MCI_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered2$ensembl_gene_id)

Other_glyco <- Other %>%
    right_join(AD_glyco[,c(1,2,4:5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
    dplyr::select(1:2, 4, 9:10, 5:8)
Other_glyco_filtered1 <- Other_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
Other_glyco_filtered2 <- Other_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered2$ensembl_gene_id)

anova_glyco <- anova %>%
    right_join(AD_glyco[,c(1,2,4:5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
    dplyr::select(1:2, 4, 8:9, 5:7)
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
addDataFrame(bm$MCI_glyco_sig_exclu[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2, row.names=FALSE)
addDataFrame(bm$Other_glyco_sig_exclu[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2+4, row.names=FALSE)
addDataFrame(bm$anova_glyco_sig_exclu[,6:8], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2+4+4, row.names=FALSE)

sheet <- createSheet(wb, "glyco_sig")
addDataFrame(bm$AD_glyco_sig, sheet=sheet, startColumn=1, row.names=T)
addDataFrame(bm$MCI_glyco_sig[, 6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig)+2, row.names=FALSE)
addDataFrame(bm$Other_glyco_sig[, 6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig)+2+4, row.names=FALSE)
addDataFrame(bm$anova_glyco_sig[,6:8], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig)+2+4+4, row.names=FALSE)

sheet <- createSheet(wb, "glyco")
addDataFrame(bm$AD_glyco, sheet=sheet, startColumn=1, row.names=T)
addDataFrame(bm$MCI_glyco[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco)+2, row.names=FALSE)
addDataFrame(bm$Other_glyco[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco)+2+4, row.names=FALSE)
addDataFrame(bm$anova_glyco[,6:8], sheet=sheet, startColumn=ncol(bm$AD_glyco)+2+4+4, row.names=FALSE)

saveWorkbook(wb, paste0("../output/rosmap_modified_paper_", region, ".xlsx"))