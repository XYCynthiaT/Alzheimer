# pkgs <- c("shiny", "dplyr", "ggplot2", "plotly", "shinydashboard", 
#           "HTSet", "DT")
# for (pkg in pkgs) {
#         suppressPackageStartupMessages(library(pkg, character.only = T))
# }

load("data/pathway.rda")

pathwayUI <- function(id) {
        ns <- NS(id)
        tagList(
                fluidRow(
                        box(
                                title = "Over Representation Test",
                                width = 6,
                                solidHeader = T,
                                status = "primary",
                                a("The over representation analysis is based on up- and down- regulated genes with adj.P.Val < 0.05."),
                                dataTableOutput(ns("enrich"))
                        ),
                        box(
                                title = "Dot Plot for Over Representation Test",
                                width = 6,
                                solidHeader = T,
                                status = "primary",
                                uiOutput(ns("control_ort")),
                                plotOutput(ns("dot_ort"))
                        )
                ),
                fluidRow(
                        box(
                                title = "Gene Set Enrichment Analysis (GSEA)",
                                width = 6,
                                solidHeader = T,
                                status = "primary",
                                dataTableOutput(ns("gse"))
                        ),
                        box(
                                title = "Dot Plot for GSEA",
                                width = 6,
                                solidHeader = T,
                                status = "primary",
                                uiOutput(ns("control_gse")),
                                plotOutput(ns("dot_gse"))
                        )
                ),
                fluidRow(
                        box(
                                title = "Visulization of GSEA Results",
                                width = 12,
                                solidHeader = T,
                                status = "primary",
                                h3(textOutput(ns("title")), align='center'),
                                plotOutput(ns("gseaPlot"))
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
                        output$control_ort <- renderUI({
                                ns <- session$ns
                                sliderInput(
                                        ns("categories"), 
                                        "The number of categories", 
                                        5, 
                                        sum(pathway[[region()]]$enrichResult$both@result$p.adjust<0.05),
                                        10, step = 1
                                )
                        })
                        output$control_gse <- renderUI({
                                ns <- session$ns
                                sliderInput(
                                        ns("categories_gse"), 
                                        "The number of categories", 
                                        5, 
                                        nrow(pathway[[region()]]$gseaResult),
                                        10, step = 1
                                )
                        })
                        output$dot_ort <- renderPlot({
                                if (region() == 'CBE') {
                                        dotplot(pathway$CBE$enrichResult$both, showCategory = input$categories)
                                } else {
                                        dotplot(pathway$TCX$enrichResult$both, showCategory = input$categories)
                                }
                        })
                        output$dot_gse <- renderPlot({
                                if (region() == 'CBE') {
                                        dotplot(pathway$CBE$gseaResult, showCategory = input$categories_gse)
                                } else {
                                        dotplot(pathway$TCX$gseaResult, showCategory = input$categories_gse)
                                }
                        })
                        output$gseaPlot <- renderPlot({
                                if (region() == 'CBE') {
                                        gseaplot2(pathway$CBE$gseaResult, 
                                                  geneSetID = input$gse_rows_selected)
                                } else {
                                        gseaplot2(pathway$TCX$gseaResult, 
                                                  geneSetID = input$gse_rows_selected)
                                }
                        })
                        output$title <- renderText({
                                if (region() == 'CBE') {
                                        return(
                                                paste0(
                                                        pathway$CBE$gseaResult@result$ID[input$gse_rows_selected],
                                                        ", ",
                                                        pathway$CBE$gseaResult@result$Description[input$gse_rows_selected]
                                                )
                                        )
                                } else {
                                        return(
                                                paste0(
                                                        pathway$TCX$gseaResult@result$ID[input$gse_rows_selected],
                                                        ", ",
                                                        pathway$TCX$gseaResult@result$Description[input$gse_rows_selected]
                                                )
                                        )
                                } 
                        })
                }
        )
}