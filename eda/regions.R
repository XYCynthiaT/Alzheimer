setwd(dirname(parent.frame(2)$ofile))

library(tidyverse)
library(readxl)
library(limma)

tcx <- read_excel("../output/mayo_TCX_glycosylation.xlsx",
                   sheet = 1, range = "B2:L130")
bm22 <- read_excel("../output/msbb_BM22.xlsx",
                   sheet = 1, range = "B2:K56")
bm36 <- read_excel("../output/msbb_BM36.xlsx",
                   sheet = 1, range = "B2:K10")
dpc <- read_excel("../output/rosmap_dorsolateral.xlsx",
                  sheet = 1, range = "B2:K41")
bm10 <- read_excel("../output/msbb_BM10.xlsx",
                   sheet = 1, range = "B2:K38")
bm44 <- read_excel("../output/msbb_BM44.xlsx",
                   sheet = 1, range = "B2:K10")
## Venn Diagram
pool1 <- c(tcx$hgnc_symbol, bm22$hgnc_symbol, bm36$hgnc_symbol, bm10$hgnc_symbol, bm44$hgnc_symbol) %>%
    unique()
venn_mat1 <- data.frame(
    tcx = pool1 %in% tcx$hgnc_symbol,
    bm22 = pool1 %in% bm22$hgnc_symbol,
    bm36 = pool1 %in% bm36$hgnc_symbol,
    bm10 = pool1 %in% bm10$hgnc_symbol,
    bm44 = pool1 %in% bm44$hgnc_symbol,
    dpc = pool1 %in% dpc$hgnc_symbol
)
vennDiagram(venn_mat1)

## tcx VS bm22&bm36
## Venn Diagram
pool <- c(tcx$hgnc_symbol, bm22$hgnc_symbol, bm36$hgnc_symbol) %>%
    unique()
venn_mat <- data.frame(
    tcx = pool %in% tcx$hgnc_symbol,
    bm22 = pool %in% bm22$hgnc_symbol,
    bm36 = pool %in% bm36$hgnc_symbol
)

vennDiagram(venn_mat)
#### tcx VS bm22
vennDiagram(venn_mat[,1:2])

tcx_bm22 <- intersect(tcx$hgnc_symbol, bm22$hgnc_symbol)
tcx2 <- filter(tcx, hgnc_symbol %in% tcx_bm22)
bm222 <- filter(bm22, hgnc_symbol %in% tcx_bm22)

## heatmap
library(pheatmap)
df <- left_join(tcx2, bm222[, c(1, 5, 7:10)], by = c("ensembl_gene_id", "pathway_id")) %>%
    select(-starts_with("pathway")) %>%
    distinct()
df_color <- data.frame(
    tcx = df$logFC.x,
    bm22 = df$logFC.y
)
rownames(df_color) <- df$hgnc_symbol
(p <- pheatmap(df_color, show_rownames = T, show_colnames = T, cluster_cols = F))
full_join(df[p$tree_row$order,c(1, 5)], bm222[,1:6], by = c("ensembl_gene_id", "hgnc_symbol"))  %>%
    View()

## BM10, 44 and Dorsolateral prefrontal cortex
pool <- c(dpc$hgnc_symbol, bm10$hgnc_symbol, bm44$hgnc_symbol) %>%
    unique()
venn_mat <- data.frame(
    dpc = pool %in% dpc$hgnc_symbol,
    bm10 = pool %in% bm10$hgnc_symbol,
    bm44 = pool %in% bm44$hgnc_symbol
)

vennDiagram(venn_mat)

dpc_bm10 <- intersect(dpc$hgnc_symbol, bm10$hgnc_symbol)

dpc2 <- filter(dpc, hgnc_symbol %in% dpc_bm10)
bm102 <- filter(bm10, hgnc_symbol %in% dpc_bm10)

# heatmap
library(pheatmap)
df <- left_join(dpc2, bm102[, c(1, 5, 7:10)], by = c("ensembl_gene_id", "pathway_id")) %>%
    select(-starts_with("pathway")) %>%
    distinct()
df_color <- data.frame(
    dpc = df$logFC.x,
    bm10 = df$logFC.y
)
rownames(df_color) <- df$hgnc_symbol
(p <- pheatmap(df_color, show_rownames = T, show_colnames = T, cluster_cols = F))
full_join(df[p$tree_row$order,1:2], dpc2[,1:6], by = c("ensembl_gene_id", "hgnc_symbol"))  %>%
    View()

## Dorsolateral prefrontal cortex and temporal cortex
pool <- c(dpc$hgnc_symbol, tcx$hgnc_symbol) %>%
    unique()
venn_mat <- data.frame(
    dpc = pool %in% dpc$hgnc_symbol,
    tcx = pool %in% tcx$hgnc_symbol
)

vennDiagram(venn_mat)

dpc_tcx <- intersect(dpc$hgnc_symbol, tcx$hgnc_symbol)

dpc3 <- filter(dpc, hgnc_symbol %in% dpc_tcx)
tcx3 <- filter(tcx, hgnc_symbol %in% dpc_tcx)

# heatmap
library(pheatmap)
df <- left_join(dpc3, tcx3[, c(1, 5, 8:11)], by = c("ensembl_gene_id", "pathway_id")) %>%
    select(-starts_with("pathway")) %>%
    distinct()
df_color <- data.frame(
    dpc = df$logFC.x,
    tcx = df$logFC.y
)
rownames(df_color) <- df$hgnc_symbol
(p <- pheatmap(df_color, show_rownames = T, show_colnames = T, cluster_cols = F))
full_join(df[p$tree_row$order,1:2], dpc3[,1:6], by = c("ensembl_gene_id", "hgnc_symbol"))  %>%
    View()
