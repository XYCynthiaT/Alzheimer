library(tidyverse)
library(org.Hs.eg.db)

glycan_o <- read.csv("../raw-data/Participating Molecules [R-HSA-5173105].tsv", sep = "\t") %>%
    mutate(Class = "O-glycosylation")
glycan_n <- read.csv("../raw-data/Participating Molecules [R-HSA-446203].tsv", sep = "\t") %>%
    mutate(Class = "N-glycosylation")

glycan <- rbind(glycan_o, glycan_n) %>%
    dplyr::select(-X) %>%
    separate(MoleculeName, into = c("id", "Name"), sep = " ", ) %>%
    dplyr::select(-id)

Ensembl_name <- mapIds(
    x = org.Hs.eg.db,
    keys = glycan$Name,
    keytype = "SYMBOL",
    column = "ENSEMBL",
    multiVals = "first"
)
Ensembl_id <- mapIds(
    x = org.Hs.eg.db,
    keys = glycan$Identifier,
    keytype = "UNIPROT",
    column = "ENSEMBL",
    multiVals = "first"
)

glycan$Ensembl <- ifelse(is.na(Ensembl_name), Ensembl_id, Ensembl_name)

glycan <- glycan[complete.cases(glycan),]

glycan$entrez <- mapIds(
    x = org.Hs.eg.db,
    keys = glycan$Ensembl,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
)

saveRDS(glycan, "glyco_reactome.rds")
