setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("dplyr","HTSet", "ggplot2", "ggpubr")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

# Plotting
logcpm <- readRDS("../data/mayo_normalized_adj.rds")
DE <- readRDS("../data/DE_mayo_adj.rds")
rownames(DE$TCX$AD) <- DE$TCX$AD$ensembl_gene_id
rownames(DE$TCX$PA) <- DE$TCX$AD$ensembl_gene_id
rownames(DE$TCX$PSP) <- DE$TCX$AD$ensembl_gene_id
rownames(DE$TCX$anova) <- DE$TCX$AD$ensembl_gene_id

genes <- c('MOG', 'GALNTL6', 'ST3GAL3', 'ST6GAL2', 'SLC35C1', 'GALNT15', 'NEU3', 'ARSG', 'B3GALNT1', 'LRP1')
plots <- vector("list", length = length(genes))
i = 1
for (s in genes) {
    ens <- featureNames(logcpm$TCX)[logcpm$TCX@fdata$hgnc_symbol == s]
    df <- data.frame(
        count = logcpm$TCX$edata[ens,],
        diagnosis = logcpm$TCX$pdata$Diagnosis
    ) 
    stat.test <- tibble(
        group1 = rep(1, 3),
        group2 = c(2, 3, 4),
        p.adj = c(DE$TCX$PA[ens,'padj'], DE$TCX$PSP[ens,'padj'], DE$TCX$AD[ens,'padj']) %>%
            signif(digits = 3)
    ) 
    plots[[i]] <- ggplot(df, aes(diagnosis, count, color = diagnosis)) +
        geom_boxplot()+
        geom_jitter(width = 0.1) +
        ylab("logcpm") +
        labs(title = paste0(s, " in TCX, ANOVA Test: p.adj=", signif(DE$TCX$anova[ens,'padj'], digits = 3))) +
        theme_bw() +
        stat_pvalue_manual(
            stat.test, 
            y.position = max(df$count), step.increase = 0.1,
            label = "p.adj"
        )
    i <- i+1
}
g=1
ens <- rownames(logcpm$TCX$edata)[logcpm$TCX$fdata$hgnc_symbol == genes[g]]
DE$TCX$AD[ens,]

## ROSMAP
logcpm <- readRDS("../data/rosmap_normalized_adj.rds")
DE <- readRDS("../data/DE_rosmap_adj.rds")
rownames(DE[[1]]$AD) <- DE[[1]]$AD$ensembl_gene_id
rownames(DE[[1]]$MCI) <- DE[[1]]$MCI$ensembl_gene_id
rownames(DE[[1]]$Other) <- DE[[1]]$Other$ensembl_gene_id
rownames(DE[[1]]$anova) <- DE[[1]]$anova$ensembl_gene_id
symbols <- c('MOG', 'GALNTL6', 'ST3GAL3', 'ST6GAL2', 'SLC35C1', 'GALNT15', 'NEU3', 'ARSG', 'B3GALNT1')
g <- 1
gene <- rownames(logcpm$edata)[logcpm$fdata$hgnc_symbol==symbols[g]]
df <- data.frame(
    count = logcpm$edata[gene,],
    diagnosis = logcpm$pdata$diagnosis
) 
stat.test <- tibble(
    group1 = rep(1, 3),
    group2 = c(2, 3, 4),
    p.adj = c(DE[[1]][["MCI"]][gene,'padj'], DE[[1]][["AD"]][gene,'padj'], DE[[1]][["Other"]][gene,'padj']) %>%
        signif(digits = 3)
) 
ggplot(df, aes(diagnosis, count, color = diagnosis)) +
    geom_boxplot()+
    geom_jitter(width = 0.1) +
    ylab("logcpm") +
    labs(title = paste0(symbols[g], " in DLPFC", ", ANOVA Test: p.adj=", signif(DE[[1]][["anova"]][gene,'padj'], digits = 3))) +
    theme_bw() +
    stat_pvalue_manual(
        stat.test, 
        y.position = max(df$count), step.increase = 0.1,
        label = "p.adj"
    )
DE[[1]][["AD"]][gene,]
## MSBB
logcpm <- readRDS("../data/msbb_normalized_adj.rds")
DE <- readRDS("../data/DE_msbb_adj.rds")
region <- "BM10"
rownames(DE[[region]]$AD) <- DE[[region]]$AD$ensembl_gene_id
rownames(DE[[region]]$transition) <- DE[[region]]$transition$ensembl_gene_id
rownames(DE[[region]]$anova) <- DE[[region]]$anova$ensembl_gene_id
symbols <- c('MOG', 'GALNTL6', 'ST6GAL2', 'ST3GAL3', 'SLC35C1', 'GALNT15', 'NEU3', 'ARSG', 'B3GALNT1')
g <- 1
gene <- rownames(logcpm[[region]]$edata)[logcpm[[region]]$fdata$hgnc_symbol == symbols[g]]
df <- data.frame(
    count = logcpm[[region]]$edata[gene,],
    diagnosis = logcpm[[region]]$pdata$diagnosis
) %>%
    filter(diagnosis != "Transition")
stat.test <- tibble(
    group1 = 1,
    group2 = c(2),
    p.adj = c(DE[[region]][["AD"]][gene,'padj']) %>%
        signif(digits = 3)
) 
ggplot(df, aes(diagnosis, count, color = diagnosis)) +
    geom_boxplot()+
    geom_jitter(width = 0.1) +
    ylab("logcpm") +
    labs(title = paste0(symbols[g], " in ", region)) +
    theme_bw() +
    stat_pvalue_manual(
        stat.test, 
        y.position = max(df$count), step.increase = 0.1,
        label = "p.adj"
    )
# DE info:
DE[[region]][["AD"]][gene,]