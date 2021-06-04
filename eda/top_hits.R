mayo_de <- readRDS("../data/DE_mayo_adj.rds")
glyco <- readRDS("../data/glyco_reactome.rds")
msbb_de <- readRDS("../data/DE_msbb_adj.rds")
rosmap_de <- readRDS("../data/DE_rosmap_adj.rds")

library(tidyverse)
library(ggplot2)
library(pheatmap)

plotVol <- function(DEtbl){
    df <- DEtbl %>%
        mutate(glyco = ifelse(ensembl_gene_id %in% glyco$Ensembl, T, F))
    ggplot(df, aes(logFC, -log(pval)))+
        geom_point(alpha=0.2)+
        geom_point(data = filter(df, glyco==T), aes(logFC, -log(pval)), color="red", alpha=0.5)+
        geom_hline(yintercept = -log(0.05))+
        theme_bw()
}

region_top_hits <- function(DEtble, top_n=50){
    DEtble %>%
        arrange(pval) %>%
        head(top_n) %>%
        .[,"ensembl_gene_id"]
}
DEtbls <- list(TCX=mayo_de$TCX$AD, CBE=mayo_de$CBE$AD, PF=msbb_de$BM10$AD, 
            STG=msbb_de$BM22$AD, PHG=msbb_de$BM36$AD, IFG=msbb_de$BM44$AD,
            DLPFC=rosmap_de$dorsolateral$AD)
sig <- lapply(DEtbls, function(x){
    filter(x, padj<0.05) 
})

sig_hits_comb <- inner_join(
    sig$TCX, sig$CBE,
    by = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description"),
    suffix=c(".TCX", ".CBE")
) %>%
    inner_join(
        sig$PF,
        by = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description")
    ) %>%
    inner_join(
        sig$STG,
        by = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description"),
        suffix=c(".PF", ".STG")
    ) %>%
    inner_join(
        sig$PHG,
        by = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description")
    ) %>%
    inner_join(
        sig$IFG, 
        by = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description"),
        suffix=c(".PHG", ".IFG")
    ) %>%
    inner_join(
        sig$DLPFC,
        by = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description")
    ) %>%
    select(c(1,2,4), starts_with("logFC"), starts_with("padj"))

intersect(sig_hits_comb$hgnc_symbol, glyco$Name) #7/197

top <- lapply(DEtbls, function(x){
    arrange(x, pval) %>%
        head(50)
})
top_hits_pos <- full_join(
    top$TCX, top$PHG,
    by = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description"),
    suffix=c(".TCX", ".PHG")
) %>%
    full_join(
        top$IFG,
        by = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description")
    ) %>%
    full_join(
        top$DLPFC,
        by = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description"),
        suffix=c(".IFG", ".DLPFC")
    ) %>%
    select(c(1,2,4), starts_with("logFC"), starts_with("padj"))
intersect(sig_hits_pos$hgnc_symbol, glyco$Name) #31/1253

top_hits <- sapply(DEtbls, region_top_hits)
top_hits2 <- as.character(top_hits)[!duplicated(as.character(top_hits))]
top_hits_detailed <- lapply(DEtbls, function(x){filter(x,ensembl_gene_id %in% top_hits2)})

top_hits_comb %>%
    select("hgnc_symbol", starts_with("logFC"), starts_with("padj")) %>%
    mutate(logFC.TCX)
    rename(logFC.DLPFC=logFC, padj.DLPFC=padj) %>%
    pivot_longer(cols = 2:ncol(.), names_to = "region", values)
plotDot <- function(count_top){
    ggplot(count_top, aes(x = region, fill = direction, y = hgnc_symbol)) +
        geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 2, binwidth = 0.5)+
        theme_bw()
}
df <- top_hits_comb %>%
    select("hgnc_symbol", starts_with("logFC")) %>%
    rename(logFC.DLPFC=logFC) %>%
    column_to_rownames("hgnc_symbol")
df <- df[complete.cases(df),]
df[df>2] <- 2
df[df< -2] <- -2
pheatmap(df, color = colorRampPalette(c("navy", "white", "red"))(50),
         show_rownames = T, show_colnames = T, 
         cluster_cols = F, na_col = "grey")

mayorna <- readRDS("../data/mayo.rds")
mayonorm <- readRDS("../data/mayo_normalized.rds")
plotBox <- function(expression, diagnosis, symbol){
    data.frame(
        Expression = expression,
        Diagnosis = diagnosis
    )%>%
        ggplot(aes(Diagnosis, Expression))+
        geom_boxplot(aes(color=Diagnosis))+
        geom_jitter(aes(color=Diagnosis), width = 0.1)+
        theme_bw()+
        labs(title = symbol, y = "logcmp")
}
