setwd(dirname(parent.frame(2)$ofile))

pkgs=c("tidyverse", "biomaRt", "HTSet")
for(pkg in pkgs){
        suppressPackageStartupMessages(library(pkg, character.only=TRUE))
}

annot = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annot = getBM(
        attributes = c(
                "ensembl_gene_id", "description", "gene_biotype", "chromosome_name",
                "start_position", "end_position", "strand"
        ),
        mart = annot
)

# CBE ---------------------------------------------------------------------

edata_cbe <- read.table(
        "../raw-data/RNAseq_geneCounts/MayoRNAseq_RNAseq_CBE_geneCounts.tsv",
        header = T
) %>%
        column_to_rownames("ensembl_id") %>%
        as.matrix()

pdata_cbe <- read.csv(
        "../raw-data/MayoRNAseq_RNAseq_CER_covariates.csv",
        header = T, 
) %>%
        filter(!is.na(Diagnosis)) %>%
        mutate(Diagnosis = sub("Pathologic Aging", "PathologicAging", Diagnosis)) %>%
        column_to_rownames("SampleID")

colnames(edata_cbe) <- sub("X", "", colnames(edata_cbe))
edata_cbe <- edata_cbe[, rownames(pdata_cbe)]

fdata_cbe <- annot %>%
        dplyr::filter(ensembl_gene_id %in% rownames(edata_cbe)) %>%
        column_to_rownames("ensembl_gene_id")
edata_cbe <- edata_cbe[rownames(fdata_cbe),]

cbe = HTSet(edata_cbe, fdata_cbe, pdata_cbe)
cbe_auto = subset_features(cbe, cbe$fdata$chromosome_name != "MT")

# TCX ---------------------------------------------------------------------

edata_tcx <- read.table(
        "../raw-data/RNAseq_geneCounts/MayoRNAseq_RNAseq_TCX_geneCounts.tsv",
        header = T
) %>%
        column_to_rownames("ensembl_id") %>%
        as.matrix()

pdata_tcx <- read.csv(
        "../raw-data/MayoRNAseq_RNAseq_TCX_covariates.csv",
        header = T, 
) %>%
        filter(!is.na(Diagnosis)) %>%
        mutate(Diagnosis = sub("Pathologic Aging", "PathologicAging", Diagnosis)) %>%
        column_to_rownames("SampleID")

colnames(edata_tcx) <- sub("X", "", colnames(edata_tcx))
edata_tcx <- edata_tcx[, rownames(pdata_tcx)]

fdata_tcx <- annot %>%
        dplyr::filter(ensembl_gene_id %in% rownames(edata_tcx)) %>%
        column_to_rownames("ensembl_gene_id")
edata_tcx <- edata_tcx[rownames(fdata_tcx),]

tcx = HTSet(edata_tcx, fdata_tcx, pdata_tcx)
tcx_auto = subset_features(tcx, tcx$fdata$chromosome_name != "MT")

RNA = list(CBE = cbe_auto, TCX = tcx_auto)
saveRDS(RNA, file = "RNA.RDS")
