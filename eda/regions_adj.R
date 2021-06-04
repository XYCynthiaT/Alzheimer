setwd(dirname(parent.frame(2)$ofile))

library(readxl)
library(VennDiagram)
library(pheatmap)
library(xlsx)
library(biomaRt)
library(tidyverse)

glyco <- readRDS("../data/glyco_reactome.rds")
glyco_names = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
glyco_names = getBM(
    attributes = c(
        "ensembl_gene_id", "hgnc_symbol", "description"
    ),
    mart = glyco_names
)
table(glyco$Class)
glyco_names <- filter(glyco_names, ensembl_gene_id %in% glyco$Ensembl)

# -------------------------------------------------------------------------

# Different in AD and Control

# -------------------------------------------------------------------------

tcx <- read_excel("../output/mayo_modified_reactome+glycolipids_TCX.xlsx",
                  sheet = 2, range = cell_cols("B:J"), col_names = T)
bm22 <- read_excel("../output/msbb_modified_reactome+glycolipids_BM22.xlsx",
                   sheet = 2, range = cell_cols("B:J"), col_names = T)
bm36 <- read_excel("../output/msbb_modified_reactome+glycolipids_BM36.xlsx",
                   sheet = 2, range = cell_cols("B:J"), col_names = T)
dpc <- read_excel("../output/rosmap_modified_reactome+glycolipids_dorsolateral.xlsx",
                  sheet = 2, range = cell_cols("B:J"), col_names = T)
bm10 <- read_excel("../output/msbb_modified_reactome+glycolipids_BM10.xlsx",
                   sheet = 2, range = cell_cols("B:J"), col_names = T)
bm44 <- read_excel("../output/msbb_modified_reactome+glycolipids_BM44.xlsx",
                   sheet = 2, range = cell_cols("B:J"), col_names = T)
cbe <- read_excel("../output/mayo_modified_reactome+glycolipids_CBE.xlsx",
                  sheet = 2, range = cell_cols("B:J"), col_names = T)

regions <- list(tcx, bm22, bm36, bm44,bm10, dpc, cbe)
names(regions) <- c('tcx', 'bm22', 'bm36', 'bm44', 'bm10', 'dpc', 'cbe')
# N-glycan
sapply(regions, function(x)nrow(filter(x, Class == "N-glycosylation")))
# O-glycan
sapply(regions, function(x)nrow(filter(x, Class == "O-glycosylation")))
# glycolipid
sapply(regions, function(x)nrow(filter(x, Class == "Glycosphingolipid")))

# -------------------------------------------------------------------------
## All regions
# -------------------------------------------------------------------------
## pool all regions together
pool_all <- c(tcx$hgnc_symbol, bm22$hgnc_symbol, 
              bm36$hgnc_symbol, bm10$hgnc_symbol, 
              bm44$hgnc_symbol, dpc$hgnc_symbol) %>%
    unique()

## glyco-genes shared in all 7 brain regions
core_all <- Reduce(intersect, list(tcx$hgnc_symbol, bm22$hgnc_symbol, 
                                   bm36$hgnc_symbol, bm10$hgnc_symbol, 
                                   bm44$hgnc_symbol, dpc$hgnc_symbol)) 
length(core_all)

## Times of apprearance of genes
count_tcx <- select(tcx, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="temporal")
count_bm22 <- select(bm22, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="bm22")
count_bm36 <- select(bm36, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="bm36")
count_dpc <- select(dpc, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="DLPFC")
count_bm10 <- select(bm10, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="bm10")
count_bm44 <- select(bm44, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="bm44")
# count_cbe <- select(cbe, c("hgnc_symbol", "logFC", "Class")) %>%
#     distinct() %>%
#     mutate(region="cerebellum")
count <- rbind(count_tcx, count_bm22, count_bm36, count_dpc, count_bm10, count_bm44)
# Split into N-/O-
count_N <- filter(count, Class=="N-glycosylation")
count_O <- filter(count, Class=="O-glycosylation")
count_L <- filter(count, Class=="Glycosphingolipid")

countSummary <- function(x){
    x %>%
        distinct() %>%
        group_by(hgnc_symbol) %>%
        summarise(present=n()) %>%
        arrange(desc(present))
}
topGenes <- function(x, count_summary, n_overlapped=4) {
    filter(x, hgnc_symbol %in% filter(count_summary, present>=n_overlapped)$hgnc_symbol) %>% 
        mutate(hgnc_symbol=factor(hgnc_symbol, levels = unique(count_summary$hgnc_symbol))) %>%
        arrange(hgnc_symbol) %>%
        mutate(region = factor(region, levels = c("DLPFC", "bm10", "bm44", "temporal", "bm22", "bm36")),
               direction = factor(ifelse(logFC>0, "up", "down"), levels = c("up", "down")))
}
# N
count_N_summary <- countSummary(count_N)
View(count_N_summary)
count_N_top <- topGenes(count_N, count_N_summary)
View(count_N_top)
filter(glyco_names, hgnc_symbol %in% count_N_top$hgnc_symbol) %>% arrange(hgnc_symbol) %>% View()
# OST genes
# OSTGenes <- c('RPN1', 'RPN2', 'DAD1', 'STT3A', 'STT3B', 'TUSC3', 'MAGT1')
# count_OST <- filter(count_N, hgnc_symbol %in% OSTGenes) %>% 
#     arrange(hgnc_symbol) %>%
#     mutate(region = factor(region, levels = c("DLPFC", "bm10", "bm44", "temporal", "bm22", "bm36")),
#            direction = factor(ifelse(logFC>0, "up", "down")))
# dot_OST <- ggplot(count_OST, aes(x = region, fill = direction, y = hgnc_symbol)) +
#     geom_dotplot(binaxis = "y", stackdir = "center", stroke=0.5, binwidth = 0.5)+
#     theme_bw()
# ggsave("dotplot_OST_ADonly.png", plot = dot_OST, device = NULL, path = "../img/",
#        scale = 1, width = 4, height = 3, units = "in",
#        dpi = 300)

# O
count_O_summary <- countSummary(count_O)
View(count_O_summary)
count_O_top <- topGenes(count_O, count_O_summary, n_overlapped = 3)
View(count_O_top)
filter(glyco_names, hgnc_symbol %in% count_O_top$hgnc_symbol) %>% arrange(hgnc_symbol) %>% View()

# L
count_L_summary <- countSummary(count_L)
View(count_L_summary)
count_L_top <- topGenes(count_L, count_L_summary, n_overlapped = 2)
View(count_L_top)
filter(glyco_names, hgnc_symbol %in% count_L_top$hgnc_symbol) %>% arrange(hgnc_symbol) %>% View()

# barplot
plotBar <- function(count_summary, annot_x = 1:6, annot_vadj = 4){
    ggplot(count_summary, aes(present))+
        geom_bar(fill="steelblue")+
        labs(x = "The number of brain regions", y = "The number of genes") +
        annotate("text", x = annot_x, y = as.numeric(table(count_summary$present))+annot_vadj, label = table(count_summary$present))+
        theme_bw()
}
# N
(bar_N <- plotBar(count_N_summary))
ggsave("../img/barplot_N_ADonly.png", plot = bar_N, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 3, units = "in",
       dpi = 300)
# O
(bar_O <- plotBar(count_O_summary))
ggsave("../img/barplot_O_ADonly.png", plot = bar_O, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 3, units = "in",
       dpi = 300)
# L
(bar_L <- plotBar(count_L_summary, annot_x = c(1:5), annot_vadj = 1))
ggsave("../img/barplot_L_ADonly.png", plot = bar_L, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 3, units = "in",
       dpi = 300)

# dotplot
plotDot <- function(count_top){
    ggplot(count_top, aes(x = region, fill = direction, y = hgnc_symbol)) +
        geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 2, binwidth = 0.5)+
        theme_bw()
}
## N-glyco
(dot_N <- plotDot(count_N_top))
ggsave("dotplot_N_ADonly.png", plot = dot_N, device = NULL, path = "../img/",
       scale = 1, width = 5, height = 6, units = "in",
       dpi = 300)
## O-glyco
(dot_O <- plotDot(count_O_top))
ggsave("dotplot_O_ADonly.png", plot = dot_O, device = NULL, path = "../img/",
       scale = 1, width = 5, height = 5, units = "in",
       dpi = 300)
## L-glyco
(dot_L <- plotDot(count_L_top))
ggsave("dotplot_L_ADonly.png", plot = dot_L, device = NULL, path = "../img/",
       scale = 1, width = 5, height = 5, units = "in",
       dpi = 300)


# -------------------------------------------------------------------------
## Temporal lobe (tcx,bm22, bm36)
# -------------------------------------------------------------------------
plotVenn <- function(class, lobe, filename){
    if (lobe=="temporal") {
        category.names = c(
            'TC', 'BM22/STG',
            'BM36/PHG'
        )
        x = list(
            filter(tcx, Class == class)$hgnc_symbol, 
            filter(bm22, Class == class)$hgnc_symbol, 
            filter(bm36, Class == class)$hgnc_symbol
        )
    } else {
        category.names = c(
            'DLPFC', 'BM10/FP',
            'BM44/IFG'
        )
        x = list(
            filter(dpc, Class == class)$hgnc_symbol, 
            filter(bm10, Class == class)$hgnc_symbol, 
            filter(bm44, Class == class)$hgnc_symbol
        )
    }
    venn.diagram(
        x = x,
        category.names = category.names,
        filename = filename,
        output = TRUE,
        imagetype = 'png',
        scaled = TRUE,
        col = 'black',
        fill = ggsci::pal_jco()(3),
        cat.col = ggsci::pal_jco()(3),
        cat.cex = 0.8,
        margrin = 0
    )
}
# N
plotVenn("N-glycosylation", "temporal", '../img/temporal_N_venn_adj_ADonly.png')
# O
plotVenn("O-glycosylation", "temporal", '../img/temporal_O_venn_adj_ADonly.png')
# L
plotVenn("Glycosphingolipid", "temporal", '../img/temporal_L_venn_adj_ADonly.png')

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
    filename = '../img/temporal_regions_venn_adj_ADonly.png',
    output = TRUE,
    imagetype = 'png',
    scaled = TRUE,
    col = 'black',
    fill = ggsci::pal_jco()(3),
    cat.col = ggsci::pal_jco()(3),
    cat.cex = 0.8,
    margrin = 0
)

core_temporal <- Reduce(intersect, list(tcx$hgnc_symbol, bm22$hgnc_symbol, 
                       bm36$hgnc_symbol) )
pool_temporal <- c(tcx$hgnc_symbol, bm22$hgnc_symbol, bm36$hgnc_symbol) %>%
    unique()
logFC_tcx <- select(tcx, c(hgnc_symbol, logFC)) %>%
    distinct() 
logFC_bm22 <- select(bm22, c(hgnc_symbol, logFC)) %>%
    distinct() 
logFC_bm36 <- select(bm36, c(hgnc_symbol, logFC)) %>%
    distinct() 

## heatmap
df_TC <- filter(logFC_tcx, hgnc_symbol %in% core_temporal) %>%
    dplyr::rename(TCX=logFC) %>%
    left_join(filter(logFC_bm22, hgnc_symbol %in% core_temporal), by="hgnc_symbol") %>%
    dplyr::rename(BM22=logFC) %>%
    left_join(filter(logFC_bm36, hgnc_symbol %in% core_temporal), by="hgnc_symbol") %>%
    dplyr::rename(BM36=logFC) %>%
    column_to_rownames("hgnc_symbol")
annot <- filter(glyco, Name %in% core_temporal) %>%
    select(Name, Class) 
df_TC_O <- df_TC[filter(annot, Class=="O-glycosylation")$Name,]
df_TC_N <- df_TC[filter(annot, Class=="N-glycosylation")$Name,]
df_TC_L <- df_TC[filter(annot, Class=="Glycosphingolipid")$Name,]

# tHE CODES BELOW ARE FOR MERGED DATA:
# duoFunGenes <- annot$Name[duplicated(annot$Name)]
# annot <- annot[!duplicated(annot$Name),]
# annot <- arrange(annot, factor(Name, levels = rownames(df_TC)))
# # p$tree_row$labels
# rownames(annot) <- annot$Name
# annot <- select(annot, 2)
# annot[duoFunGenes,] <- "both"

(p_N <- pheatmap(df_TC_N, show_rownames = T, show_colnames = T, 
               cluster_cols = F))
(p_O <- pheatmap(df_TC_O, show_rownames = T, show_colnames = T, 
                 cluster_cols = F))
(p_L <- pheatmap(df_TC_L, show_rownames = T, show_colnames = T, 
                 cluster_cols = F))

# FOR MERGED DATA:  
# temporal_annot <- full_join(filter(tcx[,c(1,2,4,5,6)], hgnc_symbol %in% core_temporal),
#                             filter(bm22[,c(1,2,5,6)], hgnc_symbol %in% core_temporal),
#                             by = c("ensembl_gene_id", "hgnc_symbol", "Class"),
#                             suffix = c(".TCX", ".BM22"))  %>%
#     full_join(filter(bm44[,c(1,2,5,6)], hgnc_symbol %in% core_temporal),
#               by = c("ensembl_gene_id", "hgnc_symbol", "Class")) %>%
#     arrange(factor(hgnc_symbol, levels = rownames(df_TC)[p$tree_row$order]))

ggsave("heatmap_temporal_N_ADonly.png", plot = p_N, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 6, units = "in",
       dpi = 300)
ggsave("heatmap_temporal_O_ADonly.png", plot = p_O, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 6, units = "in",
       dpi = 300)
ggsave("heatmap_temporal_L_ADonly.png", plot = p_L, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 6, units = "in",
       dpi = 300)

## FOR MERGED DATA:
df_TC_all <- data.frame(hgnc_symbol = pool_temporal) %>%
    left_join(logFC_tcx, by = "hgnc_symbol", ) %>%
    dplyr::rename(TCX = logFC) %>%
    left_join(logFC_bm22, by = "hgnc_symbol") %>%
    dplyr::rename(BM22 = logFC) %>%
    left_join(logFC_bm36, by = "hgnc_symbol") %>%
    dplyr::rename(BM36 = logFC) %>%
    column_to_rownames("hgnc_symbol")
annot <- filter(glyco, Name %in% pool_temporal) %>%
    select(Name, Class) 
duoFunGenes <- annot$Name[duplicated(annot$Name)]
annot <- annot[!duplicated(annot$Name),]
annot <- arrange(annot, factor(Name, levels = rownames(df_TC_all)))
rownames(annot) <- annot$Name
annot <- select(annot, 2)
annot[duoFunGenes,] <- "both"
df_TC_all[is.na(df_TC_all)] <- 0
(p_all <- pheatmap(df_TC_all, show_rownames = F, show_colnames = T,
               annotation_row = annot,
               cluster_rows = T, cluster_cols = F, na_col = "grey"))

# -------------------------------------------------------------------------
## Prefrontal lobe (BM10, 44 and Dorsolateral prefrontal cortex)
# -------------------------------------------------------------------------
## plot venn diagram
# N
plotVenn("N-glycosylation", "frontal", '../img/prefrontal_N_venn_adj_ADonly.png')
# O
plotVenn("O-glycosylation", "frontal", '../img/prefrontal_O_venn_adj_ADonly.png')
# L
plotVenn("Glycosphingolipid", "frontal", '../img/prefrontal_L_venn_adj_ADonly.png')

# MERGED
venn.diagram(
    x = list(
        dpc$hgnc_symbol, bm10$hgnc_symbol, 
        bm44$hgnc_symbol
    ),
    category.names = c(
        'DLPFC', 'BM10/anterior prefrontal cortex',
        "BM44/Broca's area"
    ),
    filename = '../img/prefrontal_regions_venn_adj_ADonly.png',
    output = TRUE,
    imagetype = 'png',
    scaled = TRUE,
    col = 'black',
    fill = ggsci::pal_jco()(3),
    cat.col = ggsci::pal_jco()(3),
    cat.cex = 0.8,
    margrin = 0
)

core_prefrontal <- Reduce(intersect, list(bm44$hgnc_symbol,bm10$hgnc_symbol, dpc$hgnc_symbol)) 
pool_prefrontal <- c(dpc$hgnc_symbol, bm10$hgnc_symbol, bm44$hgnc_symbol) %>%
    unique()
logFC_dpc <- select(dpc, c(hgnc_symbol, logFC)) %>%
    distinct() 
logFC_bm10 <- select(bm10, c(hgnc_symbol, logFC)) %>%
    distinct() 
logFC_bm44 <- select(bm44, c(hgnc_symbol, logFC)) %>%
    distinct() 

## heatmap
df_FC <- data.frame(hgnc_symbol = core_prefrontal) %>%
    left_join(logFC_dpc, by = "hgnc_symbol", ) %>%
    dplyr::rename(DLPFC = logFC) %>%
    left_join(logFC_bm10, by = "hgnc_symbol") %>%
    dplyr::rename(BM10 = logFC) %>%
    left_join(logFC_bm44, by = "hgnc_symbol") %>%
    dplyr::rename(BM44 = logFC) %>%
    column_to_rownames("hgnc_symbol")
annot <- filter(glyco, Name %in% core_prefrontal) %>%
    select(Name, Class) 
df_FC_O <- df_FC[filter(annot, Class=="O-glycosylation")$Name,]
df_FC_N <- df_FC[filter(annot, Class=="N-glycosylation")$Name,]
df_FC_L <- df_FC[filter(annot, Class=="Glycosphingolipid")$Name,]

# duoFunGenes <- annot$Name[duplicated(annot$Name)]
# annot <- annot[!duplicated(annot$Name),]
# annot <- arrange(annot, factor(Name, levels = rownames(df_FC)))
# rownames(annot) <- annot$Name
# annot <- select(annot, 2)
# annot[duoFunGenes,] <- "both"

(p_N <- pheatmap(df_FC_N, show_rownames = T, show_colnames = T, 
               cluster_cols = F))
(p_O <- pheatmap(df_FC_O, show_rownames = T, show_colnames = T, 
                 cluster_cols = F))
(p_L <- pheatmap(df_FC_L, show_rownames = T, show_colnames = T, 
                 cluster_cols = F))

# frontal_annot <- full_join(filter(dpc[,c(1,2,4,5,6)], hgnc_symbol %in% core_prefrontal), 
#                            filter(bm10[,c(1,2,5,6)], hgnc_symbol %in% core_prefrontal),
#                            by = c("ensembl_gene_id", "hgnc_symbol", "Class"),
#                            suffix = c(".DLPFC", ".BM10"))  %>%
#     full_join(filter(bm44[,c(1,2,5,6)], hgnc_symbol %in% core_prefrontal),
#               by = c("ensembl_gene_id", "hgnc_symbol", "Class"))  %>%
#     arrange(factor(hgnc_symbol, levels = rownames(df_FC)[p$tree_row$order])) 
ggsave("heatmap_frontal_N_ADonly.png", plot = p_N, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 6, units = "in",
       dpi = 300)
ggsave("heatmap_frontal_O_ADonly.png", plot = p_O, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 6, units = "in",
       dpi = 300)
ggsave("heatmap_frontal_L_ADonly.png", plot = p_L, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 6, units = "in",
       dpi = 300)


# MERGED:
df_FC_all <- data.frame(hgnc_symbol = pool_prefrontal) %>%
    left_join(logFC_dpc, by = "hgnc_symbol", ) %>%
    dplyr::rename(DLPFC = logFC) %>%
    left_join(logFC_bm10, by = "hgnc_symbol") %>%
    dplyr::rename(BM10 = logFC) %>%
    left_join(logFC_bm44, by = "hgnc_symbol") %>%
    dplyr::rename(BM44 = logFC) %>%
    column_to_rownames("hgnc_symbol")
annot <- filter(glyco, Name %in% pool_prefrontal) %>%
    select(Name, Class) 
duoFunGenes <- annot$Name[duplicated(annot$Name)]
annot <- annot[!duplicated(annot$Name),]
annot <- arrange(annot, factor(Name, levels = rownames(df_FC_all)))
rownames(annot) <- annot$Name
annot <- select(annot, 2)
annot[duoFunGenes,] <- "both"
df_FC_all[is.na(df_FC_all)] <- 0
(p_all <- pheatmap(df_FC_all, show_rownames = F, show_colnames = T,
                   annotation_row = annot,
                   cluster_rows = T, cluster_cols = F, na_col = "grey"))

# -------------------------------------------------------------------------
## Memory (Temporal cortex and Dorsolateral prefrontal cortex)
# -------------------------------------------------------------------------
## plot venn diagram
venn.diagram(
    x = list(
        dpc$hgnc_symbol, tcx$hgnc_symbol 
    ),
    category.names = c(
        'DLPFC', 'TCX'
    ),
    filename = '../img/memory_regions_venn_adj_ADonly.png',
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
core_memory <- Reduce(intersect, list(tcx$hgnc_symbol, dpc$hgnc_symbol)) 
pool_memory <- c(dpc$hgnc_symbol, tcx$hgnc_symbol) %>%
    unique()

## heatmap
df_M <- inner_join(logFC_dpc, logFC_tcx, by = "hgnc_symbol") %>%
    dplyr::rename(DLPFC = logFC.x, TCX = logFC.y) %>%
    column_to_rownames("hgnc_symbol")
annot <- filter(glyco, Name %in% core_memory) %>%
    select(Name, Class) 
duoFunGenes <- annot$Name[duplicated(annot$Name)]
annot <- annot[!duplicated(annot$Name),]
annot <- arrange(annot, factor(Name, levels = rownames(df_M)))
rownames(annot) <- annot$Name
annot <- select(annot, 2)
annot[duoFunGenes,] <- "both"

(p <- pheatmap(df_M, show_rownames = T, show_colnames = T, 
               cluster_cols = F, annotation_row = annot))
frontal_annot <- full_join(filter(dpc[,c(1,2,4,5,6)], hgnc_symbol %in% core_prefrontal), 
                           filter(bm10[,c(1,2,5,6)], hgnc_symbol %in% core_prefrontal),
                           by = c("ensembl_gene_id", "hgnc_symbol", "Class"),
                           suffix = c(".DLPFC", ".BM10"))  %>%
    full_join(filter(bm44[,c(1,2,5,6)], hgnc_symbol %in% core_prefrontal),
              by = c("ensembl_gene_id", "hgnc_symbol", "Class"))  %>%
    arrange(factor(hgnc_symbol, levels = rownames(df_FC)[p$tree_row$order])) 

memory_annot <- full_join(filter(dpc[,c(1,2,4,5,6)], hgnc_symbol %in% core_memory), 
                          filter(tcx[,c(1,2,5,6)], hgnc_symbol %in% core_memory),
                          by = c("ensembl_gene_id", "hgnc_symbol", "Class"),
                          suffix = c(".DLPFC", ".TCX"))  %>%
    arrange(factor(hgnc_symbol, levels = rownames(df_M)[p$tree_row$order])) 
ggsave("heatmap_memory_regions_adj_ADonly.png", plot = p, device = NULL, path = "../img/",
       scale = 1, width = 5, height = 7, units = "in",
       dpi = 300)
# -------------------------------------------------------------------------

# Heatmap of all regions
df_all <- inner_join(rownames_to_column(df_TC_all), rownames_to_column(df_FC_all), by = "rowname") %>%
    filter_all(all_vars(. != 0)) %>%
    column_to_rownames("rowname")
df_all_all <- full_join(rownames_to_column(df_TC_all), rownames_to_column(df_FC_all), by = "rowname") %>%
    column_to_rownames("rowname")
annot <- filter(glyco, Name %in% core_all) %>%
    select(Name, Class) 
duoFunGenes <- annot$Name[duplicated(annot$Name)]
annot <- annot[!duplicated(annot$Name),]
annot <- arrange(annot, factor(Name, levels = rownames(df_all)))
rownames(annot) <- annot$Name
annot <- select(annot, 2)
annot[duoFunGenes,] <- "both"
(p <- pheatmap(df_all, show_rownames = T, show_colnames = T, 
               cluster_cols = F, annotation_row = annot))
all_annot <- full_join(temporal_annot, frontal_annot,
                           by = c("ensembl_gene_id", "hgnc_symbol", "Class", "description")) %>%
    dplyr::rename(logFC.BM36=logFC.x, logFC.BM44=logFC.y) %>%
    arrange(factor(hgnc_symbol, levels = rownames(df_all)[p$tree_row$order])) 
ggsave("heatmap_all_regions_adj_ADonly.png", plot = p, device = NULL, path = "../img/",
       scale = 1, width = 6, height = 6, units = "in",
       dpi = 300)


wb <- createWorkbook()
sheet <- createSheet(wb, "temporal")
addDataFrame(temporal_annot, sheet=sheet, startColumn=1, row.names=T)
sheet <- createSheet(wb, "frontal")
addDataFrame(frontal_annot, sheet=sheet, startColumn=1, row.names=T)
sheet <- createSheet(wb, "memory")
addDataFrame(memory_annot, sheet=sheet, startColumn=1, row.names=T)
saveWorkbook(wb, "../img/heatmap_annot_adj.xlsx")

# -------------------------------------------------------------------------

# Exclusive to AD

# -------------------------------------------------------------------------
tcx <- read_excel("../output/mayo_modified_adj_TCX.xlsx",
                  sheet = 1, range = "B2:J94", col_names = T)
bm22 <- read_excel("../output/msbb_modified_reactome_BM22.xlsx",
                   sheet = 1, range = "B2:J56", col_names = T)
bm36 <- read_excel("../output/msbb_modified_reactome_BM36.xlsx",
                   sheet = 1, range = "B2:J20", col_names = T)
dpc <- read_excel("../output/rosmap_modified_adj_dorsolateral.xlsx",
                  sheet = 1, range = "B2:J72", col_names = T)
bm10 <- read_excel("../output/msbb_modified_reactome_BM10.xlsx",
                   sheet = 1, range = "B2:J94", col_names = T)
bm44 <- read_excel("../output/msbb_modified_reactome_BM44.xlsx",
                   sheet = 1, range = "B2:J15", col_names = T)
cbe <- read_excel("../output/mayo_modified_adj_CBE.xlsx",
                  sheet = 1, range = "B2:J49", col_names = T)

regions <- list(tcx, bm22, bm36, bm44,bm10, dpc, cbe)
names(regions) <- c('tcx', 'bm22', 'bm36', 'bm44', 'bm10', 'dpc', 'cbe')
# N
sapply(regions, function(x)nrow(filter(x, Class == "N-glycosylation")))
# O
sapply(regions, function(x)nrow(filter(x, Class == "O-glycosylation")))
## FUNCTION: select by pathways
# subset_pathways <- function(x, pathwayID = c("hsa00510", "hsa00512", "hsa00513", "hsa00514", "hsa00515")){
#     filter(x, pathway_id %in% pathwayID)
# }
# -------------------------------------------------------------------------
## All regions
# -------------------------------------------------------------------------
## pool all regions together
pool_all <- c(tcx$hgnc_symbol, bm22$hgnc_symbol, 
              bm36$hgnc_symbol, bm10$hgnc_symbol, 
              bm44$hgnc_symbol, dpc$hgnc_symbol) %>%
    unique()

## glyco-genes shared in all 7 brain regions
core_all <- Reduce(intersect, list(tcx$hgnc_symbol, bm22$hgnc_symbol, 
                                   bm36$hgnc_symbol, bm10$hgnc_symbol, 
                                   bm44$hgnc_symbol, dpc$hgnc_symbol)) 
length(core_all)
## Times of apprearance of genes
count_tcx <- select(tcx, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="temporal")
count_bm22 <- select(bm22, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="bm22")
count_bm36 <- select(bm36, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="bm36")
count_dpc <- select(dpc, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="DLPFC")
count_bm10 <- select(bm10, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="bm10")
count_bm44 <- select(bm44, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="bm44")
count_cbe <- select(cbe, c("hgnc_symbol", "logFC", "Class")) %>%
    distinct() %>%
    mutate(region="cerebellum")
count <- rbind(count_tcx, count_bm22, count_bm36, count_dpc, count_bm10, count_bm44)
# Split into N-/O-
count_N <- filter(count, Class=="N-glycosylation")
count_O <- filter(count, Class=="O-glycosylation")
# N
count_N_summary <- count_N %>%
    distinct() %>%
    group_by(hgnc_symbol) %>%
    summarise(present=n()) %>%
    arrange(desc(present))
View(count_N_summary)
# O
count_O_summary <- count_O %>%
    distinct() %>%
    group_by(hgnc_symbol) %>%
    summarise(present=n()) %>%
    arrange(desc(present))
View(count_O_summary)
# OST genes
OSTGenes <- c('RPN1', 'RPN2', 'DAD1', 'STT3A', 'STT3B', 'TUSC3', 'MAGT1')
count_OST <- filter(count_N, hgnc_symbol %in% OSTGenes) %>% 
    arrange(hgnc_symbol) %>%
    mutate(region = factor(region, levels = c("DLPFC", "bm10", "bm44", "temporal", "bm22", "bm36")),
           direction = factor(ifelse(logFC>0, "up", "down")))
dot_OST <- ggplot(count_OST, aes(x = region, fill = direction, y = hgnc_symbol)) +
    geom_dotplot(binaxis = "y", stackdir = "center", stroke=0.5, binwidth = 0.5)+
    theme_bw()
ggsave("dotplot_OST_ADonly.png", plot = dot_OST, device = NULL, path = "../img/",
       scale = 1, width = 4, height = 3, units = "in",
       dpi = 300)

# barplot
# N
bar_N <- ggplot(count_N_summary, aes(present))+
    geom_bar(fill="steelblue")+
    labs(x = "The number of brain regions", y = "The number of genes") +
    annotate("text", x = 1:4, y = as.numeric(table(count_N_summary$present))+4, label = table(count_N_summary$present))+
    theme_bw()
ggsave("barplot_N_ADspecific.png", plot = bar_N, device = NULL, path = "../img/",
       scale = 1, width = 2, height = 3, units = "in",
       dpi = 300)
# O
bar_O <- ggplot(count_O_summary, aes(present))+
    geom_bar(fill="steelblue")+
    labs(x = "The number of brain regions", y = "The number of genes") +
    annotate("text", x = 1:6, y = as.numeric(table(count_O_summary$present))+2, label = table(count_O_summary$present))+
    theme_bw()
ggsave("barplot_O_ADonly.png", plot = bar_O, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 3, units = "in",
       dpi = 300)
# dotplot
## N-glyco
count_N_top <- filter(count_N, hgnc_symbol %in% filter(count_N_summary, present>1)$hgnc_symbol) %>% 
    arrange(hgnc_symbol) %>%
    mutate(region = factor(region, levels = c("DLPFC", "bm10", "bm44", "temporal", "bm22", "bm36")),
           direction = factor(ifelse(logFC>0, "up", "down")))
dot_N <- ggplot(count_N_top, aes(x = region, fill = direction, y = hgnc_symbol)) +
    geom_dotplot(binaxis = "y", stackdir = "center", stroke=0.5, binwidth = 0.5)+
    theme_bw()
ggsave("dotplot_N_ADspec.png", plot = dot_N, device = NULL, path = "../img/",
       scale = 1, width = 5, height = 6, units = "in",
       dpi = 300)
filter(glyco_names, hgnc_symbol %in% dot_N$data$hgnc_symbol) %>%
    arrange(hgnc_symbol) %>%
    View()
## O-glyco
count_O_top <- filter(count_O, hgnc_symbol %in% filter(count_O_summary, present>1)$hgnc_symbol) %>% 
    arrange(hgnc_symbol) %>%
    mutate(region = factor(region, levels = c("DLPFC", "bm10", "bm44", "temporal", "bm22", "bm36")),
           direction = factor(ifelse(logFC>0, "up", "down")))
View(count_O_top)
dot_O <- ggplot(count_O_top, aes(x = region, fill = direction, y = hgnc_symbol)) +
    geom_dotplot(binaxis = "y", stackdir = "center", stroke=0.5, binwidth = 0.5)+
    theme_bw()
ggsave("dotplot_O_ADspec.png", plot = dot_O, device = NULL, path = "../img/",
       scale = 1, width = 5, height = 6, units = "in",
       dpi = 300)
filter(glyco_names, hgnc_symbol %in% dot_O$data$hgnc_symbol) %>%
    arrange(hgnc_symbol) %>%
    View()

# -------------------------------------------------------------------------
## Temporal lobe (tcx,bm22, bm36)
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
    filename = '../img/temporal_regions_venn_adj.png',
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
## FIX ME:
df_tcx_bm22_adj <- df_tcx_bm22
df_tcx_bm22_adj[df_tcx_bm22>0.7] <- 0.7
df_tcx_bm22_adj[df_tcx_bm22< -0.7] <- -0.7
(p <- pheatmap(df_tcx_bm22_adj, show_rownames = T, show_colnames = T, 
               cluster_cols = F))
temporal_annot <- full_join(filter(tcx[,1:7], hgnc_symbol %in% rownames(df_tcx_bm22)), 
                            filter(bm22[,c(1,2,5,7)], hgnc_symbol %in% rownames(df_tcx_bm22)),
                            by = c("ensembl_gene_id", "hgnc_symbol"),
                            suffix = c(".TCX", ".BM22"))  %>%
    arrange(factor(hgnc_symbol, levels = rownames(df_tcx_bm22)[p$tree_row$order]))

ggsave("heatmap_temporal_regions_adj.png", plot = p, device = NULL, path = "../img/",
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
    filename = '../img/prefrontal_regions_venn_adj.png',
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
                           by = c("ensembl_gene_id", "hgnc_symbol"),
                           suffix = c(".DLPFC", ".BM10"))  %>%
    arrange(factor(hgnc_symbol, levels = rownames(df_dpc_bm10)[p$tree_row$order])) 
ggsave("heatmap_frontal_regions_adj.png", plot = p, device = NULL, path = "../img/",
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
    filename = '../img/memory_regions_venn_adj.png',
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
                          by = c("ensembl_gene_id", "hgnc_symbol"),
                          suffix = c(".DLPFC", ".TCX"))  %>%
    arrange(factor(hgnc_symbol, levels = rownames(df_dpc_tcx)[p$tree_row$order])) 
ggsave("heatmap_memory_regions_adj.png", plot = p, device = NULL, path = "../img/",
       scale = 1, width = 3, height = 6, units = "in",
       dpi = 300)
wb <- createWorkbook()
sheet <- createSheet(wb, "temporal")
addDataFrame(temporal_annot, sheet=sheet, startColumn=1, row.names=T)
sheet <- createSheet(wb, "frontal")
addDataFrame(frontal_annot, sheet=sheet, startColumn=1, row.names=T)
sheet <- createSheet(wb, "memory")
addDataFrame(memory_annot, sheet=sheet, startColumn=1, row.names=T)
saveWorkbook(wb, "../img/heatmap_annot_adj.xlsx")