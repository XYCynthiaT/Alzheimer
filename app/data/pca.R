setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("dplyr", "ggplot2", "pheatmap", "ggfortify", "HTSet")
for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = T))
}
load("differentialExp.rda")
rna_normalized <- readRDS("../../data/RNA_normalized.RDS")

# volcano plots
# fit_cbe$trend <- ifelse(
#         fit_cbe$adj.P.Val>=0.05, 
#         'stable',
#         ifelse(
#                 fit_cbe$logFC >= 0,
#                 'up',
#                 'down'
#         )
# )
# fit_tcx$trend <- ifelse(
#         fit_tcx$adj.P.Val>=0.05, 
#         'stable',
#         ifelse(
#                 fit_tcx$logFC >= 0,
#                 'up',
#                 'down'
#         )
# )
# 
# vol_cbe <- ggplot(fit_cbe, aes(logFC, -log(adj.P.Val))) +
#         geom_point(aes(color = trend), alpha = 0.2) + 
#         geom_hline(yintercept = -log(0.05)) +
#         theme_bw()
# vol_tcx <- ggplot(fit_tcx, aes(logFC, -log(adj.P.Val))) +
#         geom_point(aes(color = trend), alpha = 0.2) + 
#         geom_hline(yintercept = -log(0.05)) +
#         theme_bw()


# PCA data
sig <- rownames(fit_cbe)[fit_cbe$adj.P.Val<0.05]
sig_genes_cbe <- rna_normalized$CBE[sig]
df_cbe <- sig_genes_cbe$edata %>%
        t() %>%
        as.data.frame()
df_diag_cbe = cbind(df_cbe, diagnosis = sig_genes_cbe$pdata$Diagnosis)
print("start prcomp computation for CBE.")
pca_cbe <- prcomp(df_cbe, scale. = T)
print("done with CBE.")

sig <- rownames(fit_tcx)[fit_tcx$adj.P.Val<0.05]
sig_genes_tcx <- rna_normalized$TCX[sig]
df_tcx = sig_genes_tcx$edata %>%
        t() %>%
        as.data.frame()
df_diag_tcx = cbind(df_tcx, diagnosis = sig_genes_tcx$pdata$Diagnosis)
print("start prcomp computation for TCX.")
pca_tcx <- prcomp(df_tcx, scale. = T)
print("done with TCX")

# PCA_cbe <- autoplot(prcomp(df_cbe, scale. = T), 
#                     data = df_diag_cbe, colour = 'diagnosis',
#                     frame = TRUE, frame.type = 'norm')+
#         theme_bw()


# PCA_tcx <- autoplot(prcomp(df_tcx, scale. = T), 
#                     data = df_diag_tcx, colour = 'diagnosis',
#                     frame = TRUE, frame.type = 'norm')+
#         theme_bw()

# heatmaps
# filter out genes that are insignificant
# cutoff <- rownames(fit_cbe)[fit_cbe$adj.P.Val < 0.05]
# sig <- subset_features(rna_normalized$CBE, features = cutoff)
# sig <- subset_samples(sig, samples = (sig$pdata$Diagnosis %in% c("AD", "Con")))
# df_scale <- scale(sig$edata) 
# df_scale[df_scale > 2] = 2
# df_scale[df_scale < -2] = -2
# heatcol = data.frame(group = sig$pdata$Diagnosis)
# rownames(heatcol) = colnames(df_scale)
# pheatmap(df_scale, show_rownames = F, show_colnames = F, annotation_col = heatcol)

save(pca_cbe, pca_cbe, df_diag_cbe, df_diag_tcx, file = "pca.rda")
# save(vol_cbe, vol_tcx, PCA_cbe, PCA_tcx, file = "plots.rda")
