high <- fit_cbe[abs(fit_cbe$logFC)>1000,] %>% rownames()
highRna <- rna$CBE[high,] 
highdf <- highRna$edata %>% 
        t() %>%
        cbind(diagnosis = highRna$pdata$Diagnosis) 
highdf_long <- highdf %>%
        as.data.frame() %>%
        tibble::rownames_to_column("subj") %>%
        tidyr::pivot_longer(2:(ncol(.)-1), names_to = "id", values_to = "exp")%>%
        mutate(exp = as.numeric(exp), 
               diagnosis = factor(diagnosis, levels = c("Control","PathologicAging", "PSP", "AD")))


ggplot(highdf_long, aes(diagnosis, exp, color = diagnosis)) +
        geom_boxplot() +
        geom_point(aes(group = diagnosis)) +
        theme_bw() +
        facet_wrap(~id, scales = "free")

edata_cbe <- read.table(
        "../raw-data/RNAseq_geneCounts/MayoRNAseq_RNAseq_CBE_geneCounts.tsv",
        header = T
) %>%
        tibble::column_to_rownames("ensembl_id") %>%
        as.matrix()
colnames(edata_cbe) <- sub("X", "", colnames(edata_cbe))
edata_cbe <- edata_cbe[rownames(rna$CBE$fdata), rownames(rna$CBE$pdata)]

rnaCount <- HTSet(edata_cbe, rna$CBE$fdata, rna$CBE$pdata)

MTRatio <- apply(rnaCount$edata, 2, function(x){sum(x[rnaCount$fdata$chromosome_name == "MT"])/sum(x)})
MTRatio_df <- data.frame(
        subj = names(MTRatio),
        MTRatio = MTRatio
)
ggplot(MTRatio_df) +
        geom_histogram(aes(MTRatio), bins = 40, fill = "white", color = "black") +
        theme_bw() + 
        ylab("subject counts")

tmp <- filter(MTRatio_df, MTRatio>0.25)
tmp %>% group_by(diagnosis) %>% summarise(subjCount = n())
rnaCount$pdata %>% group_by(Diagnosis) %>% summarise(count = n())

nas <- apply(highdf[2:dim(highdf)[2]], 1, function(x){sum(x == 0)})
names(nas) <- highdf$ensemblId
is.na(rna$CBE$edata["ENSG00000210082", ]) %>% sum()