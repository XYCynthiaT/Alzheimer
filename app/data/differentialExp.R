setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("limma", "dplyr","HTSet", "edgeR")
for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = T))
}

# load files from data
rna <- readRDS("RNA.RDS")

# Create a Differential Gene Expression List Object
d0_cbe <- DGEList(rna$CBE$edata)
d0_tcx <- DGEList(rna$TCX$edata)
# calculate normalization factors to scale the raw library size(number of reads)
# using the function calcNormFactors, uses TMM (weighted trimmed means of M-values
# to the reference)
d0_cbe <- calcNormFactors(d0_cbe)
d0_tcx <- calcNormFactors(d0_tcx)

# filter out genes
# hist(rowSums(cpm(d0_cbe)), xlim = c(0, 500), breaks = 100000)
cutoff_cbe <- rowSums(cpm(d0_cbe)) > 10
cutoff_tcx <- rowSums(cpm(d0_tcx)) > 10
d_cbe <- d0_cbe[cutoff_cbe,]
d_tcx <- d0_tcx[cutoff_tcx,]
# group <- factor(rna$CBE$pdata$Sex)
# plotMDS(d_cbe, col = as.numeric(group))
# plotMDS(cpm(d_cbe), col = as.numeric(group))

# extracting normalized expression table
logcpm_cbe <- cpm(d_cbe, prior.count=2, log=TRUE)
logcpm_tcx <- cpm(d_tcx, prior.count=2, log=TRUE)
fdata_cbe <- rna$CBE$fdata[rownames(logcpm_cbe),]
fdata_tcx <- rna$TCX$fdata[rownames(logcpm_tcx),]
pdata_cbe <- rna$CBE$pdata 
pdata_cbe$Diagnosis <- factor(pdata_cbe$Diagnosis, labels = c("Con", "PA", "PSP", "AD"),
                              levels = c("Control", "PathologicAging", "PSP", "AD")
)
pdata_tcx <- rna$TCX$pdata 
pdata_tcx$Diagnosis <- factor(pdata_tcx$Diagnosis, labels = c("Con", "PA", "PSP", "AD"),
                              levels = c("Control", "PathologicAging", "PSP", "AD")
)

n_cbe <- HTSet(edata = logcpm_cbe, fdata = fdata_cbe, pdata = pdata_cbe)
n_tcx <- HTSet(edata = logcpm_tcx, fdata = fdata_tcx, pdata = pdata_tcx)
n_RNA = list(CBE = n_cbe, TCX = n_tcx)

# fit the linear model
design_cbe <- model.matrix(~ 0+Diagnosis, data = rna$CBE$pdata)
transform_cbe <- voom(d_cbe, design_cbe, plot = F)
fit_cbe <- lmFit(transform_cbe, design_cbe)
contr_cbe = makeContrasts(DiagnosisAD - DiagnosisControl, levels = colnames(coef(fit_cbe)))
fit_cbe <- contrasts.fit(fit_cbe, contr_cbe) %>%
        eBayes() %>%
        topTable(sort.by = NULL, number = Inf)

design_tcx <- model.matrix(~ 0+Diagnosis, data = rna$TCX$pdata)
transform_tcx <- voom(d_tcx, design_tcx, plot = F)
fit_tcx <- lmFit(transform_tcx, design_tcx) 
contr_tcx <- makeContrasts(DiagnosisAD - DiagnosisControl, levels = colnames(coef(fit_tcx)))
fit_tcx <- contrasts.fit(fit_tcx, contr_tcx)%>%
        eBayes() %>%
        topTable(number = Inf, sort.by = NULL)

# attempt
# design_tcx2 <- model.matrix(~ Diagnosis * Gender, data = rna$TCX$pdata)
# transform_tcx2 <- voom(d_tcx, design_tcx2, plot = F)
# fit_tcx2 <- lmFit(transform_tcx2, design_tcx2) 
# fit_tcx2 <- contrasts.fit(fit_tcx2, coef = 2)%>%
#         eBayes() %>%
#         topTable(number = Inf, sort.by = NULL)
# 
# design_tcx3 <- model.matrix(~ 0 + Diagnosis + Gender, data = rna$TCX$pdata)
# transform_tcx3 <- voom(d_tcx, design_tcx3, plot = F)
# fit_tcx3 <- lmFit(transform_tcx3, design_tcx3) 
# contr_tcx3 = makeContrasts(DiagnosisAD - DiagnosisControl, levels = colnames(coef(fit_tcx3)))
# fit_tcx3 <- contrasts.fit(fit_tcx3, contr_tcx3)%>%
#         eBayes() %>%
#         topTable(number = Inf, sort.by = NULL)

# datatable
a <- fit_cbe %>%
        tibble::rownames_to_column("ensembl_id") 
b <- rna$CBE$fdata %>%
        tibble::rownames_to_column("ensembl_id") %>%
        mutate(description = sub("\\[.*\\]", "", description)) %>%
        dplyr::select(1:3)
output_lm_cbe <- left_join(a, b, by = "ensembl_id")

c <- fit_tcx %>%
        tibble::rownames_to_column("ensembl_id") 
d <- rna$TCX$fdata %>%
        tibble::rownames_to_column("ensembl_id") %>%
        mutate(description = sub("\\[.*\\]", "", description)) %>%
        dplyr::select(1:3)
output_lm_tcx <- left_join(c, d, by = "ensembl_id")

save(fit_cbe, fit_tcx, output_lm_cbe, output_lm_tcx, file = "differentialExp.rda")
saveRDS(n_RNA, "RNA_normalized.RDS")
