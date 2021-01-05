setwd(dirname(parent.frame(2)$ofile))

DE <- readRDS("../data/DE.rds")
pkgs <- c("limma", "dplyr","HTSet", "edgeR", "biomaRt", "pheatmap", "ggplot2")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

# run once----
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

targets_symbol <- c("TREM2", "APOD", "CLU", "APOE", "LRP1", "CD33", 
                    "VLDLR", "LDLR", "ABCA1", "ABCA7", "PLTP", "HMGCR", 
                    "ABCG1", "SOAT1", "SOAT2",
                    "CD36", "SCARB1", "SCARA3", "SCARA5", "MSR1",
                    "ACSS1", "ACSS2", "ACSS3", "ACSM1", "ACSM2A", "ACSM2B", 
                    "ACSM3", "ACSM4", "ACSM5", "ACSBG1", "ACSBG2", 
                    paste0("SLC27A", 1:6), paste0("ACSL", c(1, 3:6)),
                    "CPT1A", "CPT1B", "CPT1C", "CPT2", 
                    "ACADS", "ACADSB", "ACADM", "ACADL", "ACADVL",
                    "ACAD8", "ACAD9", "ACAD10", "ACAD11", "ECHS1", "ECH1", "HADH",
                    "ACAT1", "ACAT2", "NQO1", "DDIT3", "BCL2", "CASP3", "CASP8", "CASP9",
                    "HSPA5", "EIF2AK3", "ERN1", "BAX", "BAK1", "BNIP3", "PARL", 
                    "LDHA", "LDHB", "LDHC", "LDHD"
                     )
targets_dict <- getBM(
    filters="hgnc_symbol",
    attributes=c("ensembl_gene_id", "hgnc_symbol"),
    values=targets_symbol,
    mart=mart)
targets_AD <- DE$TCX$AD %>%
    filter(ensembl_gene_id %in% targets_dict$ensembl_gene_id) %>%
    left_join(targets_dict, by = "ensembl_gene_id")
targets_PA <- DE$TCX$PA %>%
    filter(ensembl_gene_id %in% targets_dict$ensembl_gene_id) %>%
    left_join(targets_dict, by = "ensembl_gene_id")
targets_PSP <- DE$TCX$PSP %>%
    filter(ensembl_gene_id %in% targets_dict$ensembl_gene_id) %>%
    left_join(targets_dict, by = "ensembl_gene_id")
targets_anova <- DE$TCX$anova %>%
    filter(ensembl_gene_id %in% targets_dict$ensembl_gene_id) %>%
    left_join(targets_dict, by = "ensembl_gene_id")
data.frame(
    ensembl_gene_id = targets_AD$ensembl_gene_id,
    hgnc_symbol = targets_AD$hgnc_symbol,
    description = targets_AD$description,
    anova.pval = targets_anova$pval,
    anova.padj = targets_anova$padj,
    AD.logFC = targets_AD$logFC,
    AD.pval = targets_AD$pval,
    AD.padj = targets_AD$padj,
    PA.logFC = targets_PA$logFC,
    PA.pval = targets_PA$pval,
    PA.padj = targets_PA$padj,
    PSP.logFC = targets_PSP$logFC,
    PSP.pval = targets_PSP$pval,
    PSP.padj = targets_PSP$padj
) %>% 
    arrange(match(hgnc_symbol, targets_symbol)) %>%
    write.csv("../output/chol_anova.csv")


library(clusterProfiler)
library(org.Hs.eg.db)
targets_dict <- bitr(targets_symbol, "SYMBOL", "ENSEMBL", org.Hs.eg.db)
colnames(targets_dict) <- c("hgnc_symbol", "ensembl_gene_id")
targets_dict <- targets_dict[!duplicated(targets_dict$ensembl_gene_id),]

# Pathway mapping
df <- data.frame(
    ensembl_gene_id = targets_AD$ensembl_gene_id,
    hgnc_symbol = targets_AD$hgnc_symbol,
    description = targets_AD$description,
    AD.logFC = targets_AD$logFC,
    AD.pval = targets_AD$pval,
    AD.padj = targets_AD$padj,
    PA.logFC = targets_PA$logFC,
    PA.pval = targets_PA$pval,
    PA.padj = targets_PA$padj,
    PSP.logFC = targets_PSP$logFC,
    PSP.pval = targets_PSP$pval,
    PSP.padj = targets_PSP$padj
) %>% 
    arrange(match(hgnc_symbol, targets_symbol))
geneList <- df$AD.logFC
names(geneList) <- id

enrich <- enrichKEGG(
    gene = names(geneList),
    organism = 'hsa'
)
library("pathview")
pathview(gene.data  = geneList,
         pathway.id = "hsa00071",
         species    = "hsa",
         limit      = list(gene=max(abs(geneList)), cpd=1))