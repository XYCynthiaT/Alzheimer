setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("dplyr", "HTSet", "xlsx")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

DE <- readRDS("../data/DE_mayo_adj.rds")
reactome <- readRDS("../data/glyco.rds") %>%
    dplyr::rename(ensembl_gene_id = Ensembl)
# reactome <- readRDS("../data/glyco_reactome.rds") %>%
#     dplyr::rename(ensembl_gene_id = Ensembl)

for (region in names(DE)) {
    AD <- DE[[region]]$AD
    PA <- DE[[region]]$PA
    PSP <- DE[[region]]$PSP
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
        filter(ensembl_gene_id %in% distinctGenes)
    # Exlude genes chaged in PSP
    sigAD <- ifelse(AD$padj>0.05, 0, ifelse(AD$logFC<0, -1, 1))
    sigPSP <- ifelse(PSP$padj>0.05, 0, ifelse(PSP$logFC<0, -1, 1))
    sigBoth <- sigAD == sigPSP
    keep <- (sigAD != 0)&!sigBoth
    # Reactome genes
    glycoGenes <- AD$ensembl_gene_id %in% reactome$ensembl_gene_id
    AD_glyco_trend <- AD %>%
        dplyr::filter(keep&glycoGenes) %>%
        left_join(distinct(reactome), by = "ensembl_gene_id") %>%
        dplyr::select(1:2, 4, 9:11, 5:8) 
    
    PA_glyco <- PA %>%
        right_join(AD_glyco[,c(1,2,4:5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
        dplyr::select(1:2, 4, 9:10, 5:8)
    PA_glyco_filtered1 <- PA_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
    PA_glyco_filtered2 <- PA_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered2$ensembl_gene_id)
    PA_glyco_trend <- PA_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_trend$ensembl_gene_id)
    
    PSP_glyco <- PSP %>%
        right_join(AD_glyco[,c(1,2,4:5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
        dplyr::select(1:2, 4, 9:10, 5:8)
    PSP_glyco_filtered1 <- PSP_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
    PSP_glyco_filtered2 <- PSP_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered2$ensembl_gene_id)
    PSP_glyco_trend <- PSP_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_trend$ensembl_gene_id)
    
    anova_glyco <- anova %>%
        right_join(AD_glyco[,c(1,2,4:5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
        dplyr::select(1:2, 4, 8:9, 5:7)
    anova_glyco_filtered1 <- anova_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
    anova_glyco_filtered2 <- anova_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered2$ensembl_gene_id)
    anova_glyco_trend <- anova_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_trend$ensembl_gene_id)
    
    bm <- list(
        AD_glyco = AD_glyco,
        AD_glyco_sig = AD_glyco_filtered1,
        AD_glyco_sig_exclu = AD_glyco_filtered2,
        AD_glyco_trend = AD_glyco_trend,
        PA_glyco = PA_glyco,
        PA_glyco_sig = PA_glyco_filtered1,
        PA_glyco_sig_exclu = PA_glyco_filtered2,
        PA_glyco_trend = PA_glyco_trend,
        PSP_glyco = PSP_glyco,
        PSP_glyco_sig = PSP_glyco_filtered1,
        PSP_glyco_sig_exclu = PSP_glyco_filtered2,
        PSP_glyco_trend = PSP_glyco_trend,
        anova_glyco = anova_glyco,
        anova_glyco_sig = anova_glyco_filtered1,
        anova_glyco_sig_exclu = anova_glyco_filtered2,
        anova_glyco_trend = anova_glyco_trend
    )
    
    wb <- createWorkbook()
    sheet <- createSheet(wb, "glyco_sig_exclu")
    addDataFrame(bm$AD_glyco_sig_exclu, sheet=sheet, startColumn=1, row.names=T)
    addDataFrame(bm$PA_glyco_sig_exclu[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2, row.names=FALSE)
    addDataFrame(bm$PSP_glyco_sig_exclu[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2+4, row.names=FALSE)
    addDataFrame(bm$anova_glyco_sig_exclu[,6:8], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2+4+4, row.names=FALSE)
    sheet <- createSheet(wb, "glyco_sig")
    addDataFrame(bm$AD_glyco_sig, sheet=sheet, startColumn=1, row.names=T)
    addDataFrame(bm$PA_glyco_sig[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig)+2, row.names=FALSE)
    addDataFrame(bm$PSP_glyco_sig[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig)+2+4, row.names=FALSE)
    addDataFrame(bm$anova_glyco_sig[,6:8], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig)+2+4+4, row.names=FALSE)
    sheet <- createSheet(wb, "glyco")
    addDataFrame(bm$AD_glyco, sheet=sheet, startColumn=1, row.names=T)
    addDataFrame(bm$PA_glyco[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2, row.names=FALSE)
    addDataFrame(bm$PSP_glyco[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2+4, row.names=FALSE)
    addDataFrame(bm$anova_glyco[,6:8], sheet=sheet, startColumn=ncol(bm$AD_glyco_sig_exclu)+2+4+4, row.names=FALSE)
    sheet <- createSheet(wb, "glyco_excToPSP")
    addDataFrame(bm$AD_glyco_trend, sheet=sheet, startColumn=1, row.names=T)
    addDataFrame(bm$PA_glyco_trend[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_trend)+2, row.names=FALSE)
    addDataFrame(bm$PSP_glyco_trend[,6:9], sheet=sheet, startColumn=ncol(bm$AD_glyco_trend)+2+4, row.names=FALSE)
    addDataFrame(bm$anova_glyco_trend[,6:8], sheet=sheet, startColumn=ncol(bm$AD_glyco_trend)+2+4+4, row.names=FALSE)
    saveWorkbook(wb, paste0("../output/mayo_modified_paper_", region, ".xlsx"))
}


