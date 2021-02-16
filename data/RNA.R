setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("tidyverse", "biomaRt", "HTSet", "data.table")
for(pkg in pkgs){
        suppressPackageStartupMessages(library(pkg, character.only=TRUE))
}

annot = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annot = getBM(
        attributes = c(
                "ensembl_gene_id", "hgnc_symbol", "description", "gene_biotype", 
                "chromosome_name",
                "start_position", "end_position", "strand"
        ),
        mart = annot
)

# Mayo --------------------------------------------------------------------

path <- "../raw-data/RNAseq_geneCounts/"
files <- list.files(path, pattern = "all_counts", full.names = T)
mayo <- vector("list", length = length(files))
for (i in seq_along(files)) {
        ## edata
        edata <- fread(files[i], header = T, sep = "\t") 
        edata_mayo <- edata[-c(1:4)] %>%
                as.data.frame() %>%
                mutate(feature = sub("\\..+", "", feature)) %>%
                column_to_rownames("feature") %>%
                as.matrix()
        ## pdata
        pdata <- read.csv(
                list.files(path, pattern = "covariates", full.names = T)[i],
                header = T, 
        ) %>%
                filter(!is.na(Diagnosis)) %>%
                mutate(Diagnosis = sub("Pathologic Aging", "PathologicAging", Diagnosis)) %>%
                mutate(
                        Diagnosis = factor(Diagnosis, labels = c("Con", "PA", "PSP", "AD"),
                                           levels = c("Control", "PathologicAging", "PSP", "AD")
                        )
                ) %>%
                column_to_rownames("SampleID")
        edata_mayo <- edata_mayo[, rownames(pdata)]
        
        ## fdata
        fdata <- annot %>%
                dplyr::filter(ensembl_gene_id %in% rownames(edata_mayo)) %>%
                filter(!duplicated(ensembl_gene_id)) %>%
                column_to_rownames("ensembl_gene_id")
        edata_mayo <- edata_mayo[rownames(fdata),]
        mayo[[i]] <-  HTSet(edata_mayo, fdata, pdata) %>%
                subset_features(.$fdata$chromosome_name != "MT")
}
# names(mayo) <- basename(files) %>% substr(6, 8)
names(mayo) <- c("CBE", "TCX")

# MSBB --------------------------------------------------------------------

path <- "../raw-data/RNAseq_MSBB/"
files <- list.files(path, pattern = "all_counts", full.names = T)
rnaseq_covariates <- read.csv(
        list.files(path, pattern = "covariates", full.names = T),
        header = T 
)
metadata <- read.csv(
        list.files(path, pattern = "metadata", full.names = T),
        header = T 
) %>%
        mutate(diagnosis = ifelse(
                CERAD<=1 & Braak<=3 & CDR<=0.5, "NCI",
                ifelse(
                        CERAD>1 & Braak>3 & CDR>0.5, "AD", "Transition"
                )
        )) %>%
        mutate(
                diagnosis = factor(diagnosis, labels = c("NCI", "Transition", "AD"))
        )

edata <- fread(files, sep = "\t", header = TRUE)
edata_msbb <- edata[-(1:4),] %>%
        as.data.frame() %>%
        mutate(feature = sub("\\..+", "", feature)) %>%
        column_to_rownames("feature") %>%
        as.matrix()
specimenID_edata <- colnames(edata_msbb) 

#### check whether there is duplication:
length(specimenID_edata)
unique(specimenID_edata) %>% length() # no

## fdata
fdata <- annot %>%
        dplyr::filter(ensembl_gene_id %in% rownames(edata_msbb)) %>%
        filter(!duplicated(ensembl_gene_id)) %>%
        column_to_rownames("ensembl_gene_id") 
edata_msbb <- edata_msbb[rownames(fdata),]

### pdata:
pdata <- rnaseq_covariates %>%
        filter(fileType == "bam") %>%
        dplyr::select(-c("synapseId", "fileName", "barcode", "fileType")) %>%
        mutate(BrodmannArea = factor(BrodmannArea)) %>%
        merge(metadata, ., by.x = "individualID", by.y = "individualIdentifier") %>%
        filter(Action == "OKay" & RIN>=4 & rRNA.rate<=0.05) %>% 
        # samples with QC actions "Remap" or "Exclude", low RIN score (<4), 
        # or relatively large rRNA rate (>5%) were removed.
        dplyr::filter(sampleIdentifier %in% specimenID_edata) %>% # 908, 118 is missing
        column_to_rownames("sampleIdentifier")
common_sample <- intersect(colnames(edata_msbb), rownames(pdata))
pdata <- pdata[common_sample,]
edata_msbb <- edata_msbb[,common_sample]

msbb <- HTSet(edata_msbb, fdata, pdata) %>%
        subset_features(.$fdata$chromosome_name != "MT")

# msbb <- vector("list", length = length(files))
# for (i in seq_along(files)) {
#         ## edata
#         edata <- read.table(
#                 files[i],
#                 header = T,
#                 sep = "\t"
#         ) %>% 
#                 column_to_rownames("Ensembl.ID") %>%
#                 as.matrix()
#         ## fdata
#         fdata <- annot %>%
#                 dplyr::filter(ensembl_gene_id %in% rownames(edata)) %>%
#                 filter(!duplicated(ensembl_gene_id)) %>%
#                 column_to_rownames("ensembl_gene_id") 
#         edata <- edata[rownames(fdata),]
#         ## pdata
#         bm <- str_split(basename(files[i]), pattern = "_", simplify = T)[1,5] %>%
#                 sub(".raw", "", .) %>%
#                 paste0("BM", .)
#         pdata <- rnaseq_covariates %>%
#                 filter(fileType == "bam") %>%
#                 filter(BrodmannArea == bm) %>%
#                 merge(metadata, ., by.x = "individualID", by.y = "individualIdentifier") %>%
#                 filter(Action == "OKay" & RIN>=4 & rRNA.rate<=0.05) %>%
#                 column_to_rownames("sampleIdentifier")
#         #samples with QC actions "Remap" or "Exclude", low RIN score (<4), 
#         # or relatively large rRNA rate (>5%) were removed.
#         common_sample <- intersect(colnames(edata), rownames(pdata))
#         pdata <- pdata[common_sample,]
#         edata <- edata[,common_sample]
#         msbb[[i]] <-  HTSet(edata, fdata, pdata) %>%
#                 subset_features(.$fdata$chromosome_name != "MT")
#         names(msbb)[i] <- bm
# }

# ROSMAP ------------------------------------------------------------------

path <- "../raw-data/RNAseq_ROSMAP/"
specimen <- read.csv(
        list.files(path, pattern = "biospecimen", full.names = T),
        header = T 
) %>%
        dplyr::select(1,2,6)
## split expression data by brain regions
files <- list.files(path, pattern = "all_counts", full.names = T)
edata <- fread(
        files,
        header = T,
        sep = "\t"
)
# edata2 <- fread(
#         files[2],
#         header = T,
#         sep = "\t"
# )

#### Check the identity of first two columns:
# identical(edata1$tracking_id, edata1$gene_id)
# identical(edata2$tracking_id, edata2$gene_id)
# identical(edata1$tracking_id, edata2$gene_id)

edata_rosmap <- edata[-(1:4),] %>%
        as.data.frame() %>%
        mutate(feature = sub("\\..+", "", feature)) %>%
        column_to_rownames("feature")
specimenID_edata <- colnames(edata_rosmap) 

#### check whether there is duplication:
length(specimenID_edata)
unique(specimenID_edata) %>% length() # no

#### find the duplicated ones:
specimenID_edata[duplicated(specimenID_edata)]
# grep("492_120515_.*", colnames(edata_rosmap))
# edata_rosmap[1:10,c(6, 459, 571)] # I will take the avg
# edata_rosmap <- rownames_to_column(edata_rosmap) %>%
#         mutate(`492_120515+m` = (`492_120515_0`+`492_120515_6`+`492_120515_7`)/3) %>%
#         dplyr::select(-c(`492_120515_0`,`492_120515_6`,`492_120515_7`)) %>%
#         column_to_rownames()
# specimenID_edata <- colnames(edata_rosmap) %>%
#         substr(1, nchar(.)-2)
# colnames(edata_rosmap) <- specimenID_edata


#### align colnames of edata to brain regions and individualID:
specimen_for_edata <- dplyr::filter(specimen, specimenID %in% specimenID_edata) %>%
        mutate(tissue = factor(tissue)) # 638, one is missing
levels(specimen_for_edata$tissue) # all are from dorsolateral prefrontal cortex
edata_rosmap <- dplyr::select(edata_rosmap, specimen_for_edata$specimenID)

#### replace the colnames of edata to individualID:
identical(specimen_for_edata$specimenID, colnames(edata_rosmap)) #false
specimen_for_edata <- arrange(specimen_for_edata, factor(specimenID, levels = colnames(edata_rosmap)))
identical(specimen_for_edata$specimenID, colnames(edata_rosmap)) # true
sum(specimen_for_edata$individualID=="") # 2 sample didn't have individualID!!!! 
colnames(edata_rosmap) <- specimen_for_edata$individualID # replace specimen ID with individual
edata_rosmap <- edata_rosmap[,colnames(edata_rosmap) != ""] %>%
        as.matrix()

## pdata_rosmap
clinical <- read.csv(
        list.files(path, pattern = "clinical.csv", full.names = T),
        header = T 
) %>%
        mutate(diagnosis = ifelse(
                pmax(cogdx, dcfdx_lv, na.rm = T) == 1, "NCI",
                ifelse(
                        pmax(cogdx, dcfdx_lv, na.rm = T)<4, "MCI",
                        ifelse(
                                pmax(cogdx, dcfdx_lv, na.rm = T)<6, "AD", "Other"
                        )
                )
        )) %>%
        mutate(
                diagnosis = factor(diagnosis, levels = c("NCI", "MCI", "AD", "Other")),
                msex = factor(msex, levels = c(1, 0), labels = c('male', 'female'))
        ) %>%
        column_to_rownames("individualID")
pdata_rosmap <- clinical[colnames(edata_rosmap),]

## fdata_rosmap
fdata_rosmap <- annot %>%
        dplyr::filter(ensembl_gene_id %in% rownames(edata_rosmap)) %>%
        filter(!duplicated(ensembl_gene_id)) %>%
        column_to_rownames("ensembl_gene_id") 
edata_rosmap <- edata_rosmap[rownames(fdata_rosmap),]
dorsolateral <- HTSet(edata_rosmap, fdata_rosmap, pdata_rosmap) %>%
        subset_features(.$fdata$chromosome_name != "MT")
rosmap <- list(dorsolateral = dorsolateral)

saveRDS(mayo, file = "mayo.rds")
saveRDS(msbb, file = "msbb.rds")
saveRDS(rosmap, file = "rosmap.rds")
