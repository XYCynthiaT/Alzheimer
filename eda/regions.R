setwd(dirname(parent.frame(2)$ofile))

library(tidyverse)
library(readxl)
library(VennDiagram)
library(pheatmap)
library(xlsx)

tcx <- read_excel("../output/mayo_modified_TCX.xlsx",
                   sheet = 1, range = "B2:K121")
bm22 <- read_excel("../output/msbb_modified_BM22.xlsx",
                   sheet = 1, range = "B2:K70")
bm36 <- read_excel("../output/msbb_modified_BM36.xlsx",
                   sheet = 1, range = "B2:K31")
dpc <- read_excel("../output/rosmap_modified_dorsolateral.xlsx",
                  sheet = 1, range = "B2:K101")
bm10 <- read_excel("../output/msbb_modified_BM10.xlsx",
                   sheet = 1, range = "B2:K127")
bm44 <- read_excel("../output/msbb_modified_BM44.xlsx",
                   sheet = 1, range = "B2:K14")
cbe <- read_excel("../output/mayo_modified_CBE.xlsx",
                  sheet = 1, range = "B2:K47")

## FUNCTION: select by pathways
subset_pathways <- function(x, pathwayID = c("hsa00510", "hsa00512", "hsa00513", "hsa00514", "hsa00515")){
    filter(x, pathway_id %in% pathwayID)
}
# -------------------------------------------------------------------------
## All regions
# -------------------------------------------------------------------------
## pool all regions together
pool_all <- c(tcx$hgnc_symbol, bm22$hgnc_symbol, 
              bm36$hgnc_symbol, bm10$hgnc_symbol, 
              bm44$hgnc_symbol, dpc$hgnc_symbol,
              cbe$hgnc_symbol) %>%
    unique()

## glyco-genes shared in all 7 brain regions
core_all <- Reduce(intersect, list(tcx$hgnc_symbol, bm22$hgnc_symbol, 
                                   bm36$hgnc_symbol, bm10$hgnc_symbol, 
                                   bm44$hgnc_symbol, dpc$hgnc_symbol,
                                   cbe$hgnc_symbol)) 
length(core_all)
## heatmap of all regions

# -------------------------------------------------------------------------
## Brain regions that have high DE glyco-genes
# -------------------------------------------------------------------------
## plot venn diagram
venn.diagram(
    x = list(
        tcx$hgnc_symbol, bm22$hgnc_symbol, 
        bm10$hgnc_symbol, dpc$hgnc_symbol
    ),
    category.names = c(
        'Temporal cortex', 'BM22/primary auditory cortex',
        'BM10/anterior prefrontal cortex',
        'DLPFC'
    ),
    filename = '../img/active_regions_venn.png',
    output = TRUE,
    imagetype = 'png',
    scaled = TRUE,
    col = 'black',
    fill = ggsci::pal_jco()(4),
    cat.col = ggsci::pal_jco()(4),
    cat.cex = 0.8,
    margrin = 0
)

# ## Display saved image
# options(repr.plot.height = 12, repr.plot.weight = 12)
# library(png)
# pp <- readPNG('../img/active_regions_venn.png')
# plot.new()
# rasterImage(pp, 0, 0, 1.1, 1.1)

core_active <- Reduce(intersect, list(tcx$hgnc_symbol, bm22$hgnc_symbol, 
                                   bm10$hgnc_symbol, dpc$hgnc_symbol)) 
length(core_active)

# -------------------------------------------------------------------------
## Brain regions on temporal lobe (tcx,bm22, bm36)
# -------------------------------------------------------------------------
## plot venn diagram
venn.diagram(
    x = list(
        tcx$hgnc_symbol, bm22$hgnc_symbol, 
        bm36$hgnc_symbol
    ),
    category.names = c(
        'Temporal cortex', 'BM22/primary auditory cortex',
        'BM36/perirhinal cortex'
    ),
    filename = '../img/temporal_regions_venn.png',
    output = TRUE,
    imagetype = 'png',
    scaled = TRUE,
    col = 'black',
    fill = ggsci::pal_jco()(3),
    cat.col = ggsci::pal_jco()(3),
    cat.cex = 0.8,
    margrin = 0
)

## Only one gene shared among all temporal lob regions
Reduce(intersect, list(tcx$hgnc_symbol, bm22$hgnc_symbol, 
                                      bm36$hgnc_symbol) )
pool_temporal <- c(tcx$hgnc_symbol, bm22$hgnc_symbol, bm36$hgnc_symbol) %>%
    unique()
logFC_tcx <- select(tcx, c(hgnc_symbol, logFC)) %>%
    distinct() 
logFC_bm22 <- select(bm22, c(hgnc_symbol, logFC)) %>%
    distinct() 
logFC_bm36 <- select(bm36, c(hgnc_symbol, logFC)) %>%
    distinct() 

## heatmap excluded BM36
df_tcx_bm22 <- inner_join(logFC_tcx, logFC_bm22, by = "hgnc_symbol") %>%
    dplyr::rename(TCX = logFC.x, BM22 = logFC.y) %>%
    column_to_rownames("hgnc_symbol")
df_tcx_bm22_adj <- df_tcx_bm22
df_tcx_bm22_adj[df_tcx_bm22>0.7] <- 0.7
df_tcx_bm22_adj[df_tcx_bm22< -0.7] <- -0.7
(p <- pheatmap(df_tcx_bm22_adj, show_rownames = T, show_colnames = T, 
               cluster_cols = F))
temporal_annot <- full_join(filter(tcx[,1:7], hgnc_symbol %in% rownames(df_tcx_bm22)), 
          filter(bm22[,c(1,2,5,7)], hgnc_symbol %in% rownames(df_tcx_bm22)),
          by = c("ensembl_gene_id", "hgnc_symbol", "pathway_id"),
          suffix = c(".TCX", ".BM22"))  %>%
    arrange(factor(hgnc_symbol, levels = rownames(df_tcx_bm22)[p$tree_row$order]))

ggsave("heatmap_temporal_regions.png", plot = p, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 6, units = "in",
       dpi = 300)
## heatmap with brain regions in temporal lobe
# df_temporal <- data.frame(hgnc_symbol = pool_temporal) %>%
#     left_join(logFC_tcx, by = "hgnc_symbol", ) %>%
#     dplyr::rename(TCX = logFC) %>%
#     left_join(logFC_bm22, by = "hgnc_symbol") %>%
#     dplyr::rename(BM22 = logFC) %>%
#     left_join(logFC_bm36, by = "hgnc_symbol") %>%
#     dplyr::rename(BM36 = logFC) %>%
#     column_to_rownames("hgnc_symbol")
# pool_temporal[c(complete.cases(df_temporal), !complete.cases(df_temporal))]
# (p <- pheatmap(df_temporal, show_rownames = F, show_colnames = T, 
#                cluster_rows = F, cluster_cols = F, na_col = "grey"))

# -------------------------------------------------------------------------
## Brain regions on prefrontal lobe (BM10, 44 and Dorsolateral prefrontal cortex)
# -------------------------------------------------------------------------
## plot venn diagram
venn.diagram(
    x = list(
        dpc$hgnc_symbol, bm10$hgnc_symbol, 
        bm44$hgnc_symbol
    ),
    category.names = c(
        'DLPFC', 'BM10/anterior prefrontal cortex',
        "BM44/Broca's area"
    ),
    filename = '../img/prefrontal_regions_venn.png',
    output = TRUE,
    imagetype = 'png',
    scaled = TRUE,
    col = 'black',
    fill = ggsci::pal_jco()(3),
    cat.col = ggsci::pal_jco()(3),
    cat.cex = 0.8,
    margrin = 0
)

## Only four gene shared among all temporal lob regions
Reduce(intersect, list(bm44$hgnc_symbol,bm10$hgnc_symbol, dpc$hgnc_symbol)) 
pool_prefrontal <- c(dpc$hgnc_symbol, bm10$hgnc_symbol, bm44$hgnc_symbol) %>%
    unique()
logFC_dpc <- select(dpc, c(hgnc_symbol, logFC)) %>%
    distinct() 
logFC_bm10 <- select(bm10, c(hgnc_symbol, logFC)) %>%
    distinct() 
logFC_bm44 <- select(bm44, c(hgnc_symbol, logFC)) %>%
    distinct() 

## heatmap excluded BM44 that has only 12 genes
df_dpc_bm10 <- inner_join(logFC_dpc, logFC_bm10, by = "hgnc_symbol") %>%
    dplyr::rename(DLPFC = logFC.x, BM10 = logFC.y) %>%
    column_to_rownames("hgnc_symbol")
df_dpc_bm10_adj <- df_dpc_bm10
df_dpc_bm10_adj[df_dpc_bm10>0.7] <- 0.7
df_dpc_bm10_adj[df_dpc_bm10< -0.7] <- -0.7
(p <- pheatmap(df_dpc_bm10_adj, show_rownames = T, show_colnames = T, 
               cluster_cols = F))
frontal_annot <- full_join(filter(dpc[,1:7], hgnc_symbol %in% rownames(df_dpc_bm10)), 
          filter(bm10[,c(1,2,5,7)], hgnc_symbol %in% rownames(df_dpc_bm10)),
          by = c("ensembl_gene_id", "hgnc_symbol", "pathway_id"),
          suffix = c(".DLPFC", ".BM10"))  %>%
    arrange(factor(hgnc_symbol, levels = rownames(df_dpc_bm10)[p$tree_row$order])) 
ggsave("heatmap_frontal_regions.png", plot = p, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 6, units = "in",
       dpi = 300)

# -------------------------------------------------------------------------
## Brain regions on memory (Temporal cortex and Dorsolateral prefrontal cortex)
# -------------------------------------------------------------------------
## plot venn diagram
venn.diagram(
    x = list(
        dpc$hgnc_symbol, tcx$hgnc_symbol 
    ),
    category.names = c(
        'DLPFC', 'TCX'
    ),
    filename = '../img/memory_regions_venn.png',
    output = TRUE,
    imagetype = 'png',
    scaled = TRUE,
    col = 'black',
    fill = ggsci::pal_jco()(2),
    cat.col = ggsci::pal_jco()(2),
    cat.cex = 0.8,
    margrin = 0
)

## 38 gene shared b/t TCX and DLPFC
Reduce(intersect, list(tcx$hgnc_symbol, dpc$hgnc_symbol)) 
pool_memory <- c(dpc$hgnc_symbol, tcx$hgnc_symbol) %>%
    unique()

## heatmap
df_dpc_tcx <- inner_join(logFC_dpc, logFC_tcx, by = "hgnc_symbol") %>%
    dplyr::rename(DLPFC = logFC.x, TCX = logFC.y) %>%
    column_to_rownames("hgnc_symbol")
df_dpc_tcx_adj <- df_dpc_tcx
df_dpc_tcx_adj[df_dpc_tcx>0.7] <- 0.7
df_dpc_tcx_adj[df_dpc_tcx< -0.7] <- -0.7
(p <- pheatmap(df_dpc_tcx_adj, show_rownames = T, show_colnames = T, 
               cluster_cols = F))
memory_annot <- full_join(filter(dpc[,1:7], hgnc_symbol %in% rownames(df_dpc_tcx)), 
          filter(tcx[,c(1,2,5,7)], hgnc_symbol %in% rownames(df_dpc_tcx)),
          by = c("ensembl_gene_id", "hgnc_symbol", "pathway_id"),
          suffix = c(".DLPFC", ".TCX"))  %>%
    arrange(factor(hgnc_symbol, levels = rownames(df_dpc_tcx)[p$tree_row$order])) 
ggsave("heatmap_memory_regions.png", plot = p, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 6, units = "in",
       dpi = 300)
wb <- createWorkbook()
sheet <- createSheet(wb, "temporal")
addDataFrame(temporal_annot, sheet=sheet, startColumn=1, row.names=T)
sheet <- createSheet(wb, "frontal")
addDataFrame(frontal_annot, sheet=sheet, startColumn=1, row.names=T)
sheet <- createSheet(wb, "memory")
addDataFrame(memory_annot, sheet=sheet, startColumn=1, row.names=T)
saveWorkbook(wb, "../img/heatmap_annot.xlsx")