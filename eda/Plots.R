setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("dplyr","HTSet", "ggplot2", "ggpubr")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

# Plotting
logcpm <- readRDS("../data/tcx_normalized.rds")
rna <- readRDS("../data/mayo.rds")
DE <- readRDS("../data/DE.rds")
rownames(DE$TCX$AD) <- DE$TCX$AD$ensembl_gene_id
rownames(DE$TCX$PA) <- DE$TCX$AD$ensembl_gene_id
rownames(DE$TCX$PSP) <- DE$TCX$AD$ensembl_gene_id
rownames(DE$TCX$anova) <- DE$TCX$AD$ensembl_gene_id

# gene list
genes <- c('ENSG00000121578', 'ENSG00000134910',
           'ENSG00000253710', 'ENSG00000141429', 
           'ENSG00000101346', 'ENSG00000121578',
           'ENSG00000144057', 'ENSG00000126091', 
           'ENSG00000131386',
           'ENSG00000100979', 'ENSG00000130164',
           'ENSG00000165029', 'ENSG00000147852',
           'ENSG00000123384', 'ENSG00000120885',
           'ENSG00000057252')
plots <- vector("list", length = length(genes))
i = 1
for (ens in genes) {
    df <- data.frame(
        count = rna$TCX$edata[ens,],
        diagnosis = rna$TCX$pdata$Diagnosis
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
        labs(title = paste0(ens, ", ANOVA Test: p.adj=", signif(DE$TCX$anova[ens,'padj'], digits = 3))) +
        theme_bw() +
        stat_pvalue_manual(
            stat.test, 
            y.position = max(df$count), step.increase = 0.1,
            label = "p.adj"
        )
    i <- i+1
}

# ST3GAL3, GALNT11, GALNT17
## Mayo
genes <- c('ENSG00000126091', 
           'ENSG00000178234',
           'ENSG00000185274')
df <- data.frame(
    count = logcpm[genes[1],],
    diagnosis = rna$TCX$pdata$Diagnosis
) 
stat.test <- tibble(
    group1 = rep(1, 3),
    group2 = c(2, 3, 4),
    p.adj = c(DE$TCX$PA[genes[1],'padj'], DE$TCX$PSP[genes[1],'padj'], DE$TCX$AD[genes[1],'padj']) %>%
        signif(digits = 3)
) 
ST3GAL3 <- ggplot(df, aes(diagnosis, count, color = diagnosis)) +
    geom_boxplot()+
    geom_jitter(width = 0.1) +
    ylab("logcpm") +
    labs(title = paste0("ST3GAL3 in TCX", ", ANOVA Test: p.adj=", signif(DE$TCX$anova[genes[1],'padj'], digits = 3))) +
    theme_bw() +
    stat_pvalue_manual(
        stat.test, 
        y.position = max(df$count), step.increase = 0.1,
        label = "p.adj"
    )
## ROSMAP
logcpm <- readRDS("../data/rosmap_normalized.rds")
rna <- readRDS("../data/rosmap.rds")
DE <- readRDS("../data/DE_rosmap.rds")
rownames(DE[[1]]$AD) <- DE[[1]]$AD$ensembl_gene_id
rownames(DE[[1]]$MCI) <- DE[[1]]$MCI$ensembl_gene_id
rownames(DE[[1]]$Other) <- DE[[1]]$Other$ensembl_gene_id
rownames(DE[[1]]$anova) <- DE[[1]]$anova$ensembl_gene_id
df <- data.frame(
    count = logcpm[genes[1],],
    diagnosis = rna[[1]]$pdata$diagnosis
) 
stat.test <- tibble(
    group1 = rep(1, 3),
    group2 = c(2, 3, 4),
    p.adj = c(DE[[1]][["MCI"]][genes[1],'padj'], DE[[1]][["AD"]][genes[1],'padj'], DE[[1]][["Other"]][genes[1],'padj']) %>%
        signif(digits = 3)
) 
ST3GAL3 <- ggplot(df, aes(diagnosis, count, color = diagnosis)) +
    geom_boxplot()+
    geom_jitter(width = 0.1) +
    ylab("logcpm") +
    labs(title = paste0("ST3GAL3 in DLPFC", ", ANOVA Test: p.adj=", signif(DE[[1]][["anova"]][genes[1],'padj'], digits = 3))) +
    theme_bw() +
    stat_pvalue_manual(
        stat.test, 
        y.position = max(df$count), step.increase = 0.1,
        label = "p.adj"
    )
## MSBB
logcpm <- readRDS("../data/msbb_normalized.rds")
rna <- readRDS("../data/msbb.rds")
DE <- readRDS("../data/DE_msbb.rds")
region <- "BM44"
rownames(DE[[region]]$AD) <- DE[[region]]$AD$ensembl_gene_id
rownames(DE[[region]]$transition) <- DE[[region]]$transition$ensembl_gene_id
rownames(DE[[region]]$anova) <- DE[[region]]$anova$ensembl_gene_id
df <- data.frame(
    count = logcpm[[region]][genes[2],],
    diagnosis = rna[[region]]$pdata$diagnosis
) 
stat.test <- tibble(
    group1 = rep(1, 2),
    group2 = c(2, 3),
    p.adj = c(DE[[region]][["transition"]][genes[1],'padj'], DE[[region]][["AD"]][genes[1],'padj']) %>%
        signif(digits = 3)
) 
ggplot(df, aes(diagnosis, count, color = diagnosis)) +
    geom_boxplot()+
    geom_jitter(width = 0.1) +
    ylab("logcpm") +
    labs(title = paste0("ST3GAL3 in Broca's area", ", ANOVA Test: p.adj=", signif(DE[[region]][["anova"]][genes[1],'padj'], digits = 3))) +
    theme_bw() +
    stat_pvalue_manual(
        stat.test, 
        y.position = max(df$count), step.increase = 0.1,
        label = "p.adj"
    )
