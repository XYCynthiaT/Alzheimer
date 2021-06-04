setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("limma", "dplyr","HTSet", "edgeR", "biomaRt")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

# load files from data
rna <- readRDS("../data/msbb.rds")
# run once----
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#-------------


# -------------------------------------------------------------------------
# Split dataset by brain regions
# -------------------------------------------------------------------------
rna <- list(
    BM10 = subset_samples(rna, rna$pdata$BrodmannArea == "BM10"),
    BM22 = subset_samples(rna, rna$pdata$BrodmannArea == "BM22"),
    BM36 = subset_samples(rna, rna$pdata$BrodmannArea == "BM36"),
    BM44 = subset_samples(rna, rna$pdata$BrodmannArea == "BM44")
)

# edgeR -------------------------------------------------------------------

# region <- 'BM44'
DE_msbb <- vector("list", length(rna))
logcpm_msbb <- vector("list", length(rna))
names(DE_msbb) <- names(rna)
names(logcpm_msbb) <- names(rna)
for (region in names(rna)) {
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
    design <- model.matrix(~diagnosis+sex, data = rna[[region]]$pdata)
    colnames(design) <- c("NCI", "Transition", "AD", "male")
    # Estimate the dispersion
    d <- normalization(region, geneType = "protein_coding")
    d <- estimateDisp(d, design, robust=TRUE)
    #-----
    # d$common.dispersion
    # plotBCV(d)
    #------
    # extracting normalized expression table
    logcpm <- cpm(d, prior.count=0.5, log=TRUE)
    logcpm_pf <- rna[[region]][rownames(logcpm),]
    logcpm_pf$edata <- logcpm
    logcpm_msbb[[region]] <- logcpm_pf
    # Fitting NB model
    fit <- glmQLFit(d, design, robust=TRUE)
    # plotQLDisp(fit)
    qlfAnova <- glmQLFTest(fit, coef = 2:3)
    qlfTransition <- glmQLFTest(fit, coef = 2)
    qlfAD <- glmQLFTest(fit, coef = 3)

    # DE Tables
    resAnova <- topTags(qlfAnova, n = Inf, sort.by = "none")
    resAD <- topTags(qlfAD, n = Inf, sort.by = "none")
    resTransition <- topTags(qlfTransition, n = Inf, sort.by = "none")

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
    transition <- finalTable(region, resTransition)

    anov <- cbind(AD[,1:4], resAnova$table[,4:6]) %>%
        dplyr::rename(stat = F, pval = PValue, padj =FDR)
    rownames(anov) <- NULL

    # Genes different among groups
    deGenes <- anov$ensembl_gene_id[anov$padj<0.05]

    findDistinct <- function(){
        # return genes distinct to AD
        vennTransition <- decideTests(qlfTransition)
        vennAD <- decideTests(qlfAD)
        venn <- cbind(vennTransition, vennAD)
        # vennDiagram(venn)
        distinctAD <- (venn[, 2] != 0)&(venn[,1] != venn[,2])
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
        transition = transition
    )
    DE_msbb[[region]] <- bm
}


saveRDS(DE_msbb, "../data/DE_msbb_adj.rds")
saveRDS(logcpm_msbb, "../data/msbb_normalized_adj.rds")
