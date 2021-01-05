setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("limma", "dplyr","HTSet", "edgeR", "biomaRt")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

DE <- readRDS("data/DE.rds")
nodes <- read.csv("raw-data/iRefIndex default node.csv")
node_symbol <- nodes$hgnc
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
convt <- getBM(
    filters="hgnc_symbol",
    attributes=c("ensembl_gene_id", "hgnc_symbol"),
    values=node_symbol,
    mart=mart)

nodes_AD <- DE$TCX$AD %>%
    filter(ensembl_gene_id %in% convt$ensembl_gene_id) %>%
    left_join(convt, by = "ensembl_gene_id")
nodes_expr <- left_join(nodes, nodes_AD, by = c("hgnc" = "hgnc_symbol"))
write.csv(nodes_expr, "AD_expr.csv")
