pkgs <- c("shiny", "dplyr", "ggplot2", "plotly", "shinydashboard", "ggfortify")
for (pkg in pkgs) {
        suppressPackageStartupMessages(library(pkg, character.only = T))
}

load("data/plots.rda")

pcaUI <- function(id) {
        ns <- NS(id)
        tagList(
                column(
                        6,
                        box(
                                title = "Volcano Plot",
                                width = NULL,
                                solidHeader = T,
                                status = "primary",
                                plotlyOutput(ns("volcano"))
                        )
                ),
                column(
                        6, 
                        box(
                                title = "PCA Plot",
                                width = NULL,
                                solidHeader = T,
                                status = "primary",
                                plotlyOutput(ns("PCA"))
                        )
                )
        )
}

pcaServer <- function(id, region) {
        moduleServer(
                id,
                function(input, output, session) {
                        output$volcano <- renderPlotly({
                                if (region() == "CBE") {
                                        vol_cbe
                                } 
                                else {
                                        vol_tcx
                                }
                        })
                        output$PCA <- renderPlotly(
                                if (region() == "CBE") {
                                        PCA_cbe
                                } else {
                                        PCA_tcx
                                }
                        )
                }
        )
}