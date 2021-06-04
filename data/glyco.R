library(tidyverse)
library(org.Hs.eg.db)
library(readxl)
library(biomaRt)

# From Reactome
glycan_o <- read.csv("../raw-data/Participating Molecules [R-HSA-5173105].tsv", sep = "\t") %>%
    mutate(Class = "O-glycosylation")
glycan_n <- read.csv("../raw-data/Participating Molecules [R-HSA-446203].tsv", sep = "\t") %>%
    mutate(Class = "N-glycosylation")
glycan_l0 <- read.csv("../raw-data/Participating Molecules [R-HSA-428157].tsv", sep = "\t") %>%
    mutate(Class = "Glycosphingolipid")
glycan_l <- glycan_l0 %>% 
    separate(Molecule.Name, into = c("Molecule", "Name"), sep = " ") %>%
    dplyr::select(!c("X", "Molecule"))

# combine N-glycosylation, O-glycosylation and glycolipid
glycan <- rbind(glycan_o, glycan_n) %>%
  dplyr::select(-X) %>%
  separate(MoleculeName, into = c("id", "Name"), sep = " ") %>%
  dplyr::select(-id) %>%
  rbind(glycan_l)

# mapping the ensembl_id by uniProt_id
Ensembl_id <- mapIds(
  x = org.Hs.eg.db,
  keys = glycan$Identifier,
  keytype = "UNIPROT",
  column = "ENSEMBL",
  multiVals = "first"
)

# mapping the ensembl_id by hgnc_symbol
Ensembl_name <- mapIds(
  x = org.Hs.eg.db,
  keys = glycan$Name,
  keytype = "SYMBOL",
  column = "ENSEMBL",
  multiVals = "first"
)

# adding the ensembl_id
glycan$Ensembl <- ifelse(is.na(Ensembl_name), Ensembl_id, Ensembl_name)

glycan <- glycan[complete.cases(glycan),]

# adding entrez_id
glycan$entrez <- mapIds(
    x = org.Hs.eg.db,
    keys = glycan$Ensembl,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
)

# add description
glyco_names <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
glyco_names <- getBM(
  attributes = c(
    "ensembl_gene_id", "hgnc_symbol", "description"
  ),
  mart = glyco_names
)
filter(glyco_names, ensembl_gene_id %in% glycan$Ensembl) %>% dim()
filter(glyco_names, hgnc_symbol %in% glycan$Name) %>% dim()
description <- filter(glyco_names, hgnc_symbol %in% glycan$Name) %>%
  dplyr::select(c("hgnc_symbol", "description")) %>%
  distinct() %>%
  dplyr::rename(Name=hgnc_symbol) %>%
  mutate(description=sub("\\[.*\\]", "", description))
glycan <- left_join(glycan, description, by="Name")

# saveRDS(glycan, "glyco_reactome.rds")

# From Review Paper
glycosyltransferase <- read_excel("../raw-data/Global view of human protein glycosylation pathways and functions_suppl.xlsx", 
                                  sheet = 1, skip = 1) %>%
  mutate(
    Specificity = c(rep("Specific", 120), rep("Non-specific", 53), rep("Specific", 40)),
    Class = c(
      rep("GPI-anchor", 6),
      rep("Glycolipid", 8),
      rep("N-Glycan", 28),
      rep("O-GalNAc", 26),
      rep("O-Fucose", 6),
      rep("O-Mannose POMT-directed", 12),
      rep("O-Mannose TMTC-directed", 4),
      rep("C-Mannose", 4),
      rep("O-Glucose", 6),
      rep("O-Xylose", 16),
      rep("Hxl-Galactose", 2),
      rep("O-GlcNAc", 2),
      rep("Elongation.branching", 18),
      rep("Capping", 35),
      rep("Sulfotransferases", 40)
    )
  )
# mapping ensembl_id by uniProt_id
Ensembl_id2 <- mapIds(
  x = org.Hs.eg.db,
  keys = glycosyltransferase$`UniProt Accession`,
  keytype = "UNIPROT",
  column = "ENSEMBL",
  multiVals = "first"
)

# Gene labels use HGNC names, 
# diverging only for the GALNT family due to errors with the HGNC names - 
# GALNT20 (HGNC: GALNTL5), GALNT19 (HGNC:GALNT17), GALNT17 (HGNC: GALNTL6), 
# and are not presented in italics to ease readibility of labels.
hgnc_symbol <- mapIds(
  x = org.Hs.eg.db,
  keys = glycosyltransferase$`UniProt Accession`,
  keytype = "UNIPROT",
  column = "SYMBOL",
  multiVals = "first"
)

# glycan_enzyme
# check:
identical(names(Ensembl_id2), glycosyltransferase$`UniProt Accession`)
# if yes, bind directly, if not, use join()
glycosyltransferase <- glycosyltransferase %>%
  mutate(Ensembl=Ensembl_id2) %>%
  dplyr::select(Specificity, Class, !Pathway, everything())%>%
  dplyr::select(1:5, 11) %>%
  dplyr::rename(Name = "Gene (HGNC)", description = "Protein name (UniProt)", Identifier = "UniProt Accession")

# glycan_nonenzyme
# remove glycosyltransferases from the reactome list
non_glycosyltransferase <- filter(glycan, !(Identifier %in% glycosyltransferase$Identifier)) %>%
  mutate(Specificity = NA) %>%
  dplyr::select(8, 4, 3, 7, 5, 2)

# rbind
glyco <- rbind(glycosyltransferase, non_glycosyltransferase)

# save
saveRDS(glyco, "glyco.rds")
