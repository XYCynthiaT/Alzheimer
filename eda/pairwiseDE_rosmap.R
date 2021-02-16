setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("limma", "dplyr","HTSet", "edgeR")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

# load files from data
rna <- readRDS("../data/rosmap.rds")
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
design <- model.matrix(~diagnosis, data = rna[[region]]$pdata)
colnames(design) <- c("NCI", "MCI", "AD", "Other")
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
saveRDS(logcpm, "../data/rosmap_normalized.rds")

# Fitting NB model
fit <- glmQLFit(d, design, robust=TRUE)
# plotQLDisp(fit)
qlfAnova <- glmQLFTest(fit, coef = 2:4)
qlfMCI <- glmQLFTest(fit, coef = 2)
qlfAD <- glmQLFTest(fit, coef = 3)
qlfOther <- glmQLFTest(fit, coef = 4)

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
    rename(stat = F, pval = PValue, padj =FDR)
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
DE_rosmap <- list(dorsolateral = bm)

saveRDS(DE_rosmap, "../data/DE_rosmap.rds")
