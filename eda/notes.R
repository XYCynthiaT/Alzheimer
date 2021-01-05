setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("limma", "dplyr","HTSet", "edgeR")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

# load files from data
rna <- readRDS("RNA.RDS")

design_cbe <- model.matrix(~ 0+Diagnosis, data = rna$CBE$pdata)
fit_cbe <- lmFit(rna$CBE$edata, design_cbe)

hist(fit_cbe$Amean)
plotSA(fit_cbe)
keep <- fit_cbe$Amean > 100
fit2 <- eBayes(fit_cbe[keep,], trend=TRUE)
plotSA(fit2)


x <- rep(LETTERS[1:4], each = 10)
y <- c(rnorm(10), rnorm(10, 0.2, 1), rnorm(10, 0, 0.2), rnorm(10, 5, 2))
matr <- matrix(y, nrow = 1)
x <- factor(x)
design <- model.matrix(~0+x)
colnames(design) <- LETTERS[1:4]
fit <- lmFit(matr, design)
contrast.matrix <- makeContrasts(D-A, C-A, B-A, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix) %>%
    eBayes()
topTable(fit2, coef = 1)
topTable(fit2, coef = 2)
topTable(fit2, coef = 3)
results <- decideTests(fit2)
vennDiagram(results)

design2 <- model.matrix(~x)
colnames(design2) <- LETTERS[1:4]
fit3 <- lmFit(matr, design2) %>%
    eBayes()
topTable(fit3, coef = 2)
topTable(fit3, coef = 3)
topTable(fit3, coef = 4)
results3 <- decideTests(fit3)
vennDiagram(results3)
