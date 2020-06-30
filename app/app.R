pkgs <- c("shiny", "limma", "dplyr", "ggplot2", "plotly", "shinydashboard", "Metabase", "edgeR")
for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = T))
}

# load files from data
rna <- readRDS("../data/RNA.RDS")

# discard genes detected in less than 10% of the subjects(<16)
discard_cbe <- apply(rna$CBE$conc_table, 1, function(x) sum(x != 0))<16
discard_tcx <- apply(rna$TCX$conc_table, 1, function(x) sum(x != 0))<16

rna$CBE <- subset_features(rna$CBE, features = !discard_cbe)
rna$TCX <- subset_features(rna$TCX, features = !discard_tcx)


# Create a Differential Gene Expression List Object
d_cbe <- DGEList(rna$CBE$conc_table)
d_tcx <- DGEList(rna$TCX$conc_table)
# claculate normalization factors to scale the raw library size(number of reads)
# using the function calcNormFactors, uses TMM (weighted trimmed means of M-values
# to the reference)
d_cbe <- calcNormFactors(d_cbe)
d_tcx <- calcNormFactors(d_tcx)

# fit the linear model
# design_cbe <- model.matrix(~ Diagnosis, data = as(rna$CBE$sample_table, "data.frame"))
# fit_cbe <- lmFit(rna$CBE$conc_table, design_cbe) %>% 
#         eBayes() %>%
#         topTable(coef = "DiagnosisControl", number = Inf, sort.by = "P")
# 
# design_tcx <- model.matrix(~ Diagnosis, data = as(rna$TCX$sample_table, "data.frame"))
# fit_tcx <- lmFit(rna$TCX$conc_table, design_tcx) %>% 
#         eBayes() %>%
#         topTable(coef = "DiagnosisControl", number = Inf, sort.by = "P")

design_cbe = model.matrix(~0 + Diagnosis, data = as(rna$CBE$sample_table, "data.frame"))
transform_cbe = voom(d_cbe, design_cbe, plot = F)

fit_cbe = lmFit(transform_cbe, design_cbe)
contr = makeContrasts(DiagnosisAD - DiagnosisControl, levels = colnames(coef(fit_cbe)))
fit_cbe = contrasts.fit(fit_cbe, contr) %>%
        eBayes() %>%
        topTable(sort.by = "P", n = Inf)

design_tcx = model.matrix(~0 + Diagnosis, data = as(rna$TCX$sample_table, "data.frame"))
transform_tcx = voom(d_tcx, design_tcx, plot = F)

fit_tcx = lmFit(transform_tcx, design_tcx)
contr_tcx = makeContrasts(DiagnosisAD - DiagnosisControl, levels = colnames(coef(fit_tcx)))
fit_tcx = contrasts.fit(fit_tcx, contr_tcx) %>%
        eBayes() %>%
        topTable(sort.by = "P", n = Inf)

fit_cbe$trend <- ifelse(
        fit_cbe$adj.P.Val>=0.05, 
        'stable',
        ifelse(
                fit_cbe$logFC >= 0.1,
                'up',
                ifelse(
                        fit_cbe$logFC <= -0.1,
                        'down',
                        'stable'
                )
        )
)
fit_tcx$trend <- ifelse(
        fit_tcx$adj.P.Val>=0.05, 
        'stable',
        ifelse(
                fit_tcx$logFC >= 0.1,
                'up',
                ifelse(
                        fit_tcx$logFC <= -0.1,
                        'down',
                        'stable'
                )
        )
)

vol_cbe <- ggplot(fit_cbe, aes(logFC, -log(adj.P.Val))) +
        geom_point(aes(color = trend), alpha = 0.2) + 
        geom_hline(yintercept = -log(0.05)) +
        theme_bw()
vol_tcx <- ggplot(fit_tcx, aes(logFC, -log(adj.P.Val))) +
        geom_point(aes(color = trend), alpha = 0.2) + 
        geom_hline(yintercept = -log(0.05)) +
        theme_bw()


# PCA
df_cbe = cpm(d_cbe$counts, log = T) %>%
        t() %>%
        as.data.frame()
df_diag_cbe = bind_cols(df_cbe, diagnosis = rna$CBE$sample_table$Diagnosis)

PCA_cbe <- autoplot(prcomp(df_cbe, scale. = T), 
         data = df_diag_cbe, colour = 'diagnosis',
         frame = TRUE, frame.type = 'norm')+
        theme_bw()

df_tcx = cpm(d_tcx$counts, log = T) %>%
        t() %>%
        as.data.frame()
df_diag_tcx = bind_cols(df_tcx, diagnosis = rna$TCX$sample_table$Diagnosis)

PCA_tcx <- autoplot(prcomp(df_tcx, scale. = T), 
                    data = df_diag_tcx, colour = 'diagnosis',
                    frame = TRUE, frame.type = 'norm')+
        theme_bw()

# UI
ui <- dashboardPage(
        dashboardHeader(title = "MayoRNAseq"),
        dashboardSidebar(
                sidebarMenu(
                        menuItem("Gene expression", tabName = "gene"),
                        menuItem("Genetic variants", tabName = "variants"),
                        menuItem("Proteomics", tabName = "proteomics")
                )
        ),
        dashboardBody(
                tabItems(
                        #first tab content
                        tabItem(
                                tabName = "gene",
                                column(3, 
                                       radioButtons("region", "Brain Regions", 
                                                   choices = c("Cerebellum" = "cbe", 
                                                               "Temporal cortex" = "tcx"),
                                                   selected = "cbe"
                                       )
                                ),
                                column(9,
                                       tags$h3("Volcano Plot"),
                                       plotlyOutput('volcano'),
                                       tags$h3("PCA Plot"),
                                       plotlyOutput("PCA")
                                )
                        ),
                        #second tab content
                        tabItem(
                                tabName = "variants"
                        ),
                        #third tab content
                        tabItem(
                                tabName = "proteomics"
                        )
                )
                
        )
)

# server
server <- function(input, output){
        output$volcano <- renderPlotly({
                if (input$region == 'cbe') {
                        vol_cbe
                } else {
                        vol_tcx
                }
        })
        output$PCA <- renderPlotly({
                if (input$region == 'cbe') {
                        PCA_cbe
                } else {
                        PCA_tcx
                }
        })
}

shinyApp(ui, server)