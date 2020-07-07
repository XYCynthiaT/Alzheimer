setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("dplyr", "clusterProfiler", "biomaRt", "ggplot2")
for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = T))
}
load("data.rda")

# Cerebellum --------------------------------------------------------------

# get entrezid
# genes1 <- bitr(rownames(fit_cbe), fromType = "ENSEMBL",
#                    toType = c("ENTREZID", "SYMBOL"),
#                    OrgDb = org.Hs.eg.db) # making mistakes when matching 
#                                          # ensembl id to entrez id
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes1 <- getBM(
        filters="ensembl_gene_id",
        attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
        values=rownames(fit_cbe),
        mart=mart)

# check multiple mapping
# duplicatedID <- genes1[duplicated(genes1$ensembl_gene_id),]$ensembl_gene_id
# duplicatedEntry <- filter(genes1, ensembl_gene_id %in% duplicatedID) %>% dplyr::arrange(ensembl_gene_id)

# discard multiple entrez id, NA
genes1 <- genes1[!duplicated(genes1$ensembl_gene_id),]
genes1 <- genes1[complete.cases(genes1),]

# KEGG over-representation test (KEGG Enrichment Analysis of a gene set).
# up & down
# adj.P.val < 0.05
gene <- rownames(fit_cbe)[fit_cbe$adj.P.Val < 0.05]
gene_entrezid_cbe <- filter(genes1, ensembl_gene_id %in% gene)$entrezgene_id
kk_cbe <- enrichKEGG(gene      = gene_entrezid_cbe,
                  organism     = 'hsa')
# all_cbe_both <- kk_cbe@result[, 2:7]
glycan_cbe_both <- kk_cbe@result[grep("glycan", kk_cbe@result$Description, ignore.case = TRUE),2:7]
# view the results: arrange(all_cbe_both, pvalue) %>% head()
# sum(kk_cbe@result$p.adjust < 0.05)

# up
# logFC > 0, adj.P.val < 0.05
geneup <- rownames(fit_cbe)[fit_cbe$logFC > 0 & fit_cbe$adj.P.Val < 0.05]
geneup_entrezid_cbe <- filter(genes1, ensembl_gene_id %in% geneup)$entrezgene_id
kkup_cbe <- enrichKEGG(gene       = geneup_entrezid_cbe,
                     organism     = 'hsa')
glycan_cbe_up <- kkup_cbe@result[grep("glycan", kkup_cbe@result$Description, ignore.case = TRUE),2:7] 

# down
# logFC < 0, adj.Pvale < 0.05
genedown <- rownames(fit_cbe)[fit_cbe$logFC < 0 & fit_cbe$adj.P.Val < 0.05]
genedown_entrezid_cbe <- filter(genes1, ensembl_gene_id %in% genedown)$entrezgene_id
kkdown_cbe <- enrichKEGG(gene     = genedown_entrezid_cbe,
                     organism     = 'hsa')
glycan_cbe_down <- kkdown_cbe@result[grep("glycan", kkdown_cbe@result$Description, ignore.case = TRUE),2:7] 

# KEGG Gene Set Enrichment Analysis (Gene Set Enrichment Analysis of KEGG)
tmp <- fit_cbe %>% 
        tibble::rownames_to_column("ensembl_gene_id") %>%
        dplyr::select(1:2)
preGeneList <- merge(tmp, genes1, by = "ensembl_gene_id")
geneList_cbe <- preGeneList$logFC
names(geneList_cbe) <- preGeneList$entrezgene_id
geneList_cbe <- base::sort(geneList_cbe, decreasing = T)

kk2_cbe <- gseKEGG(
        geneList = geneList_cbe,
        organism = 'hsa',
        nPerm = 1000,
        minGSSize = 120,
        # pvalueCutoff = 0.05,
        verbose = FALSE
)
# view the results: kk2_cbe@result[, 1:10] %>% arrange(pvalue) %>% head()
# gseaplot(kk2_cbe, geneSetID = "hsa05205")

# Temporal Cortex ---------------------------------------------------------

# get entrezid
genes2 <- getBM(
        filters="ensembl_gene_id",
        attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
        values=rownames(fit_tcx),
        mart=mart)
# discard multiple entrez id, NA
genes2 <- genes2[!duplicated(genes2$ensembl_gene_id),]
genes2 <- genes2[complete.cases(genes2),]

# up & down
# adj.P.val < 0.05
gene <- rownames(fit_tcx)[fit_tcx$adj.P.Val < 0.05]
gene_entrezid_tcx <- filter(genes2, ensembl_gene_id %in% gene)$entrezgene_id
kk_tcx <- enrichKEGG(gene      = gene_entrezid_tcx,
                     organism     = 'hsa')
# all_tcx_both <- kk_tcx@result[, 2:7]
glycan_tcx_both <- kk_tcx@result[grep("glycan", kk_tcx@result$Description, ignore.case = TRUE),2:7]

# up
# logFC > 0, adj.P.val < 0.05
geneup <- rownames(fit_tcx)[fit_tcx$logFC > 0 & fit_tcx$adj.P.Val < 0.05]
geneup_entrezid_tcx <- filter(genes2, ensembl_gene_id %in% geneup)$entrezgene_id
kkup_tcx <- enrichKEGG(gene       = geneup_entrezid_tcx,
                       organism     = 'hsa')
glycan_tcx_up <- kkup_tcx@result[grep("glycan", kkup_tcx@result$Description, ignore.case = TRUE),2:7] 

# down
# logFC < 0, adj.Pvale < 0.05
genedown <- rownames(fit_tcx)[fit_tcx$logFC < 0 & fit_tcx$adj.P.Val < 0.05]
genedown_entrezid_tcx <- filter(genes2, ensembl_gene_id %in% genedown)$entrezgene_id
kkdown_tcx <- enrichKEGG(gene     = genedown_entrezid_tcx,
                         organism     = 'hsa')
glycan_tcx_down <- kkdown_tcx@result[grep("glycan", kkdown_tcx@result$Description, ignore.case = TRUE),2:7] 

# KEGG Gene Set Enrichment Analysis (Gene Set Enrichment Analysis of KEGG)
tmp <- fit_tcx %>% 
        tibble::rownames_to_column("ensembl_gene_id") %>%
        dplyr::select(1:2)
preGeneList <- merge(tmp, genes2, by = "ensembl_gene_id")
geneList_tcx <- preGeneList$logFC
names(geneList_tcx) <- preGeneList$entrezgene_id
geneList_tcx <- base::sort(geneList_tcx, decreasing = T)

kk2_tcx <- gseKEGG(
        geneList = geneList_tcx,
        organism = 'hsa',
        nPerm = 1000,
        minGSSize = 120,
        # pvalueCutoff = 0.05,
        verbose = FALSE
)

pathway <- list(
        CBE = list(
                enrichResult = list(
                        both = kk_cbe,
                        up = kkup_cbe,
                        down = kkdown_cbe
                ),
                gseaResult = kk2_cbe,
                glycan = list(
                        both = glycan_cbe_both,
                        up = glycan_cbe_up,
                        down = glycan_cbe_down
                )
        ),
        TCX = list(
                enrichResult = list(
                        both = kk_tcx,
                        up = kkup_tcx,
                        down = kkdown_tcx
                ),
                gseaResult = kk2_tcx,
                glycan = list(
                        both = glycan_tcx_both,
                        up = glycan_tcx_up,
                        down = glycan_tcx_down
                )
        )
)
save(pathway, file = "pathway.rda")

# dotplot_glycan = function(ke){
#         ke@result[grep("glycan", ke@result$Description, ignore.case = TRUE),] %>%
#                 mutate(
#                         GeneRatio = sapply(GeneRatio, function(x) eval(parse(text = x)))
#                 ) %>%
#                 arrange(GeneRatio) %>%
#                 mutate(Description = factor(Description, levels = Description)) %>%
#                 ggplot() +
#                 geom_point(aes(x = GeneRatio, y = Description, color = pvalue, size = Count)) +
#                 scale_color_gradient(low = pal_lancet()(2)[1], high = pal_lancet()(2)[2]) +
#                 labs(y = NULL) +
#                 theme_bw() +
#                 theme(
#                         axis.text = element_text(color = "black",size = 12),
#                         plot.title = element_text(hjust = 0.5)
#                 )
# }
# 
# 
# 
# # barplot
# graphics::barplot(kk_cbe, showCategory = 20)
# 
# # gene-concept network
# kk_cbeSymbol <- setReadable(kk_cbe, 'org.Hs.eg.db', 'ENTREZID')
# cnetplot(kk_cbeSymbol, foldChange=geneList_cbe)
# 
# # enrichment map
# emapplot(kk_cbe)
# 
# # browseKEGG, Alzheimer disease
# browseKEGG(kk, 'hsa05010')
# 
# # pathview
# library("pathview")
# pathview(gene.data  = geneList_cbe,
#                      pathway.id = "hsa00190",
#                      species    = "hsa",
#                      limit      = list(gene=max(abs(geneList_cbe)), cpd=1))