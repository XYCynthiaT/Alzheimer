pkgs <- c("shiny", "limma", "dplyr", "ggplot2", "plotly", "shinydashboard", 
          "HTSet", "edgeR", "DT", "ggfortify")
for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = T))
}

rna <- readRDS("../data/RNA_normalized.RDS")
load("data/differentialExp.rda")
load("data/pathway.rda")
load("data/plots.rda")

# UI
sidebar <- dashboardSidebar(
        sidebarMenu(
                # first menu
                menuItem("Gene expression", tabName = "gene",
                         menuSubItem("Differential Expression", tabName = "rna-de"),
                         menuSubItem("Volcano and PCA", tabName = "rna-pca"),
                         menuSubItem("Pathway Enrichment", tabName = "rna-pathway")
                ),
                # second menu
                menuItem("Genetic variants", tabName = "variants"),
                # third menu
                menuItem("Proteomics", tabName = "proteomics"),
                radioButtons("region", "Brain Regions",
                             choices = c("Cerebellum" = "cbe",
                                         "Temporal cortex" = "tcx"),
                             selected = "cbe"
                )
        )
)
body <- dashboardBody(
        tags$link(href="styles.css", rel="stylesheet"),
        tabItems(
                # the first subtab content
                tabItem(
                        tabName = "rna-de",
                        column(
                                6,
                                box(
                                        title = "Differential Expression",
                                        width = NULL,
                                        solidHeader = T,
                                        status = "primary",
                                        dataTableOutput("de")
                                )
                        ),
                        column(
                                6,
                                box(
                                        title = "Boxplot",
                                        width = NULL,
                                        solidHeader = T,
                                        status = "primary",
                                        collapsible = T,
                                        tags$a("Note: If the boxplot only presents one of the two regions, it means the gene expression is less than 10 count per million in the other region."),
                                        plotlyOutput("boxplot")
                                )
                        )
                ),
                # the second subtab content
                tabItem(
                        tabName = "rna-pca",
                        column(
                                6,
                                box(
                                        title = "Volcano Plot",
                                        width = NULL,
                                        solidHeader = T,
                                        status = "primary",
                                        plotlyOutput("volcano")
                                )
                        ),
                        column(
                                6, 
                                box(
                                        title = "PCA Plot",
                                        width = NULL,
                                        solidHeader = T,
                                        status = "primary",
                                        plotlyOutput("PCA")
                                )
                        )
                ),
                # the third subtab content
                tabItem(
                        tabName = "rna-pathway",
                        fluidRow(
                                box(
                                        title = "Over Representation Analysis",
                                        width = 12,
                                        solidHeader = T,
                                        status = "primary",
                                        h3("The over representation analysis is based on up- and down- regulated genes with adj.P.Val < 0.05."),
                                        dataTableOutput("kk1")
                                )
                        ),
                        fluidRow(
                                box(
                                        title = "Gene Set Enrichment Analysis",
                                        width = 12,
                                        solidHeader = T,
                                        status = "primary",
                                        dataTableOutput("kk2")
                                )
                        )
                ),
                # the forth tab content
                tabItem(tabName = "variants"),
                # the fifth tab content
                tabItem(tabName = "proteomics")
        )
)
ui <- dashboardPage(
        dashboardHeader(title = "MayoRNAseq"),
        sidebar, body
)

# server
server <- function(input, output){
        output$de <- renderDataTable({
                if (input$region == 'cbe') {
                        datatable(output_lm_cbe, 
                                  rownames = FALSE,
                                  selection = list(mode = 'single', selected = 1)
                        ) %>%
                                formatSignif(columns = 2:7, digits = 2)
                } else {
                        datatable(output_lm_tcx, 
                                  rownames = FALSE,
                                  selection = list(mode = 'single', selected = 1)
                        ) %>%
                                formatSignif(columns = 2:7, digits = 2)
                }
        })
        output$boxplot <- renderPlotly({
                if (input$region == 'cbe') {
                        id <- as.character(output_lm_cbe$ensembl_id[input$de_rows_selected])
                        if ((id %in% rownames(rna$TCX$edata))) {
                                gene_cbe <- rna$CBE[id]
                                gene_tcx <- rna$TCX[id]
                                df_cbe <- data.frame(
                                        expression = gene_cbe$edata[1,],
                                        diagnosis = gene_cbe$pdata$Diagnosis,
                                        tissue = gene_cbe$pdata$Tissue
                                )
                                df_tcx <- data.frame(
                                        expression = gene_tcx$edata[1,],
                                        diagnosis = gene_tcx$pdata$Diagnosis, 
                                        tissue = gene_tcx$pdata$Tissue
                                )
                                df <- rbind(df_cbe, df_tcx)
                                ggplot(df, aes(diagnosis, expression, color = diagnosis)) +
                                        # stat_boxplot(geom = "errorbar") +
                                        geom_boxplot() +
                                        geom_jitter(alpha = 0.3) +
                                        labs(title = id) +
                                        facet_wrap(~tissue) +
                                        theme_bw()
                        } else {
                                gene <- rna$CBE[id]
                                df <- data.frame(
                                        expression = gene$edata[1,],
                                        diagnosis = gene$pdata$Diagnosis
                                )
                                ggplot(df, aes(diagnosis, expression, color = diagnosis)) +
                                        # stat_boxplot(geom = "errorbar") +
                                        geom_boxplot() +
                                        geom_jitter(alpha = 0.3) +
                                        labs(title = paste0(id, ", ", input$region)
                                        ) +
                                        theme_bw()
                        }
                } else {
                        id <- as.character(output_lm_tcx$ensembl_id[input$de_rows_selected])
                        if ((id %in% rownames(rna$CBE$edata))) {
                                gene_cbe <- rna$CBE[id]
                                gene_tcx <- rna$TCX[id]
                                df_cbe <- data.frame(
                                        expression = gene_cbe$edata[1,],
                                        diagnosis = gene_cbe$pdata$Diagnosis,
                                        tissue = gene_cbe$pdata$Tissue
                                )
                                df_tcx <- data.frame(
                                        expression = gene_tcx$edata[1,],
                                        diagnosis = gene_tcx$pdata$Diagnosis, 
                                        tissue = gene_tcx$pdata$Tissue
                                )
                                df <- rbind(df_cbe, df_tcx)
                                ggplot(df, aes(diagnosis, expression, color = diagnosis)) +
                                        # stat_boxplot(geom = "errorbar") +
                                        geom_boxplot() +
                                        geom_jitter(alpha = 0.3) +
                                        labs(title = id) +
                                        facet_wrap(~tissue) +
                                        theme_bw()
                        } else {
                                gene <- rna$TCX[id]
                                df <- data.frame(
                                        expression = gene$edata[1,],
                                        diagnosis = gene$pdata$Diagnosis
                                )
                                ggplot(df, aes(diagnosis, expression, color = diagnosis)) +
                                        # stat_boxplot(geom = "errorbar") +
                                        geom_boxplot() +
                                        geom_jitter(alpha = 0.3) +
                                        labs(title = paste0(id, ", ", input$region)
                                        ) +
                                        theme_bw()
                        }
                }
        })
        output$volcano <- renderPlotly({
                if (input$region == 'cbe') {
                        h3("Cerebellum")
                        vol_cbe
                } else {
                        h3("Temperal Cortex")
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
        output$kk1 <- renderDataTable({
                if (input$region == 'cbe') {
                        datatable(
                                pathway$CBE$enrichResult$both@result[,2:7] %>%
                                        tibble::rownames_to_column("ID"),
                                rownames = F,
                                selection = list(mode = 'single', selected = 1)
                        ) %>%
                                formatSignif(columns = 5:7, digits = 2)
                } else {
                        datatable(
                                pathway$TCX$enrichResult$both@result[,2:7] %>%
                                        tibble::rownames_to_column("ID"),
                                rownames = F,
                                selection = list(mode = 'single', selected = 1)
                        ) %>%
                                formatSignif(columns = 5:7, digits = 2)
                }
        })
        output$kk2 <- renderDataTable({
                if (input$region == 'cbe') {
                        datatable(
                                pathway$CBE$gseaResult@result[, 1:10],
                                rownames = F,
                                selection = list(mode = 'single', selected = 1)
                        ) %>%
                                formatSignif(columns = 4:8, digits = 2)
                } else {
                        datatable(
                                pathway$TCX$gseaResult@result[, 1:10],
                                rownames = F,
                                selection = list(mode = 'single', selected = 1)
                        ) %>%
                                formatSignif(columns = 4:8, digits = 2)
                }
        })
}

shinyApp(ui, server)