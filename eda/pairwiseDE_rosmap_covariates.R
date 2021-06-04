setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("limma", "dplyr","HTSet", "edgeR", "xlsx")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

# load files from data
rna <- readRDS("../data/rosmap.rds")
rna$dorsolateral$pdata$apoe4 <- ifelse(rna$dorsolateral$pdata$apoe_genotype %in% c(24, 34, 44), "present", "absent")
# run once----
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#-------------

# edgeR -------------------------------------------------------------------

region <- "dorsolateral"
normalization <- function(region, geneType = NULL){
    group <- rna[[region]]$pdata$diagnosis
    if (!is.null(geneType)){
        dataset <- subset_features(rna[[region]], 
                                   rna[[region]]$fdata$gene_biotype %in% geneType)
        d0 <- DGEList(counts = dataset$edata, group = group)
    } else {
        d0 <- DGEList(counts = rna[[region]]$edata, group = group)
    }
    keep <- filterByExpr(d0) 
    d0 <- d0[keep, , keep.lib.sizes=FALSE]
    calcNormFactors(d0)
}
# group <- factor(rna$CBE$pdata$Sex)
# plotMDS(d_cbe, col = as.numeric(group))

# Design the model
design <- model.matrix(~diagnosis+msex+diagnosis*apoe4, data = rna[[region]]$pdata)
# colnames(design) <- c("NCI", "MCI", "AD", "Other", "female")
# Estimate the dispersion
d <- normalization(region, geneType = "protein_coding")
d <- estimateDisp(d, design, robust=TRUE)
#-----
# d$common.dispersion
# plotBCV(d)
# plotMeanVar(d, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)
#------
# # extracting normalized expression table
logcpm <- cpm(d, prior.count=0.5, log=TRUE)

# Fitting NB model
fit <- glmQLFit(d, design, robust=TRUE)
# plotQLDisp(fit)
qlfAnova <- glmQLFTest(fit, coef = 7:9)
qlfMCI <- glmQLFTest(fit, coef = 7)
qlfAD <- glmQLFTest(fit, coef = 8)
qlfOther <- glmQLFTest(fit, coef = 9)

# DE Tables
resAnova <- topTags(qlfAnova, n = Inf, sort.by = "none")
resAD <- topTags(qlfAD, n = Inf, sort.by = "none")
resMCI <- topTags(qlfMCI, n = Inf, sort.by = "none")
resOther <- topTags(qlfOther, n = Inf, sort.by = "none")

standardTable <- function(table){
    # table = table from topTags
    data.frame(
        ensembl_gene_id = rownames(table),
        logFC = table$logFC,
        stat = table$F,
        pval = table$PValue,
        padj = table$FDR
    )
}
finalTable <- function(region = region, res = resAD){
    stdTable <- res$table %>%
        standardTable()
    fdata <- rna[[region]]$fdata[stdTable$ensembl_gene_id,]
    data.frame(
        ensembl_gene_id = stdTable$ensembl_gene_id,
        hgnc_symbol = fdata$hgnc_symbol,
        gene_biotype = fdata$gene_biotype,
        description = fdata$description
    ) %>%
        bind_cols(stdTable[,-1]) 
}

AD <- finalTable(region, resAD)
MCI <- finalTable(region, resMCI)
Other <- finalTable(region, resOther)

anov <- cbind(AD[,1:4], resAnova$table[,5:7]) %>%
    dplyr::rename(stat = F, pval = PValue, padj =FDR)
rownames(anov) <- NULL

# Genes different among groups
deGenes <- anov$ensembl_gene_id[anov$padj<0.05]

findDistinct <- function(){
    # return genes distinct to AD
    vennMCI <- decideTests(qlfMCI)
    vennOther <- decideTests(qlfOther)
    vennAD <- decideTests(qlfAD)
    venn <- cbind(vennMCI, vennOther, vennAD)
    # vennDiagram(venn)
    distinctAD <- (venn[, 3] != 0)&(venn[,1] != venn[,3])&(venn[,2] != venn[,3])
    # venn[distinctAD, ] %>%
    #     as.matrix() %>%
    #     pheatmap(show_rownames = F, show_colnames = T, cluster_rows = F, cluster_cols = F)
    rownames(venn)[distinctAD]
}
distinctGenes <- findDistinct() %>% intersect(deGenes)

bm <- list(
    distinctGenes = distinctGenes,
    anova = anov,
    AD = AD,
    MCI = MCI,
    Other = Other
)
DE <- list(dorsolateral = bm)

reactome <- readRDS("../data/glyco_reactome.rds") %>%
    dplyr::rename(ensembl_gene_id = Ensembl)

region <- "dorsolateral"
AD <- DE[[region]]$AD
MCI <- DE[[region]]$MCI
Other <- DE[[region]]$Other
anova <- DE[[region]]$anova

# DE of all glycosylation related genes
AD_glyco <- AD %>%
    dplyr::filter(ensembl_gene_id %in% reactome$ensembl_gene_id) %>%
    right_join(distinct(reactome[,4:5]), by = "ensembl_gene_id") %>%
    dplyr::select(1:4, 9, 5:8)
# DE of glycosylation related genes with padj < 0.05
AD_glyco_filtered1 <- AD_glyco %>%
    dplyr::filter(padj < 0.05) 
# DE of glycosylation related genes with padj < 0.05 and distinct to AD
distinctGenes <- DE[[region]]$distinctGenes
AD_glyco_filtered2 <- AD_glyco %>%
    dplyr::filter(ensembl_gene_id %in% distinctGenes)

MCI_glyco <- MCI %>%
    right_join(AD_glyco[,c(1,2,5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
    dplyr::select(1:4, 9, 5:8)
MCI_glyco_filtered1 <- MCI_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
MCI_glyco_filtered2 <- MCI_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered2$ensembl_gene_id)

Other_glyco <- Other %>%
    right_join(AD_glyco[,c(1,2,5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
    dplyr::select(1:4, 9, 5:8)
Other_glyco_filtered1 <- Other_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered1$ensembl_gene_id)
Other_glyco_filtered2 <- Other_glyco %>%
    dplyr::filter(ensembl_gene_id %in% AD_glyco_filtered2$ensembl_gene_id)

anova_glyco <- anova %>%
    right_join(AD_glyco[,c(1,2,5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
    dplyr::select(1:4, 8, 5:7)
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

saveWorkbook(wb, paste0("../output/rosmap_modified_apoe4_", region, ".xlsx"))

wb <- createWorkbook()
sheet <- createSheet(wb, "AD")
addDataFrame(DE$dorsolateral$AD, sheet=sheet, startColumn=1, row.names=T)
sheet <- createSheet(wb, "MCI")
addDataFrame(DE$dorsolateral$MCI, sheet=sheet, startColumn=1, row.names=T)
sheet <- createSheet(wb, "Other")
addDataFrame(DE$dorsolateral$Other, sheet=sheet, startColumn=1, row.names=T)
sheet <- createSheet(wb, "ANOVA")
addDataFrame(DE$dorsolateral$anova, sheet=sheet, startColumn=1, row.names=T)
saveWorkbook(wb, paste0("../output/rosmap_modified_apoe4_whole_genes", region, ".xlsx"))


# saveRDS(DE_rosmap, "../data/DE_rosmap_adj.rds")
# saveRDS(logcpm, "../data/rosmap_normalized_adj.rds")