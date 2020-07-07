pkgs <- c("shiny", "limma", "dplyr", "ggplot2", "plotly", "shinydashboard", 
          "HTSet", "DT")
for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = T))
}

load("data/pathway.rda")

pathwayUI <- function(id) {
        ns <- NS(id)
        tagList(
                fluidRow(
                        box(
                                title = "Over Representation Analysis",
                                width = 12,
                                solidHeader = T,
                                status = "primary",
                                h3("The over representation analysis is based on up- and down- regulated genes with adj.P.Val < 0.05."),
                                dataTableOutput(ns("enrich"))
                        )
                ),
                fluidRow(
                        box(
                                title = "Enrichment Analysis",
                                width = 12,
                                solidHeader = T,
                                status = "primary",
                                dataTableOutput(ns("gse"))
                        )
                )
        )
}

pathwayServer <- function(id, region) {
        moduleServer(
                id, 
                function(input, output, session) {
                        output$enrich <- renderDataTable({
                                if (region() == 'CBE') {
                                        datatable(
                                                pathway$CBE$enrichResult$both@result[,2:7] %>%
                                                        tibble::rownames_to_column("ID"),
                                                rownames = F,
                                                selection = list(mode = 'single', selected = 1)
                                        ) %>%
                                                formatSignif(columns = 5:7, digits = 3)
                                } else {
                                        datatable(
                                                pathway$TCX$enrichResult$both@result[,2:7] %>%
                                                        tibble::rownames_to_column("ID"),
                                                rownames = F,
                                                selection = list(mode = 'single', selected = 1)
                                        ) %>%
                                                formatSignif(columns = 5:7, digits = 3)
                                }
                        })
                        output$gse <- renderDataTable({
                                if (region() == 'CBE') {
                                        datatable(
                                                pathway$CBE$gseaResult@result[,1:10],
                                                rownames = F,
                                                selection = list(mode = 'single', selected = 1)
                                        ) %>%
                                                formatSignif(columns = 4:8, digits = 3)
                                } else {
                                        datatable(
                                                pathway$TCX$gseaResult@result[,1:10],
                                                rownames = F,
                                                selection = list(mode = 'single', selected = 1)
                                        ) %>%
                                                formatSignif(columns = 4:8, digits = 3)
                                }
                        })
                }
        )
}