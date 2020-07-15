# pkgs <- c("shiny", "dplyr", "ggplot2", "plotly", "shinydashboard", "ggfortify")
# for (pkg in pkgs) {
#         suppressPackageStartupMessages(library(pkg, character.only = T))
# }

# pcaEnv = new.env(parent = .GlobalEnv)
load("data/pca.rda")
load("data/differentialExp.rda")

volcanoPlot <- function(data){
        # data are fit_xxx from differential expression analysis with linear model.
        data$trend <- ifelse(
                data$adj.P.Val>=0.05, 
                'stable',
                ifelse(
                        data$logFC >= 0,
                        'up',
                        'down'
                )
        )
        ggplot(data, aes(logFC, -log(adj.P.Val))) +
        geom_point(aes(color = trend), alpha = 0.2) + 
        geom_hline(yintercept = -log(0.05)) +
        theme_bw()
}

pcaPlot <- function(pcaData, data) {
        autoplot(pcaData, 
                data = data, colour = 'diagnosis',
                frame = TRUE, frame.type = 'norm')+
        theme_bw()
}

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
                                # if(is.null(pcaEnv$vol_cbe)){
                                #         print("loading plots.rda from volcano")
                                #         load("data/plots.rda", envir = pcaEnv)
                                #         print("done reading plots.rda from volcano")
                                # }
                                if (region() == "CBE") {
                                        volcanoPlot(fit_cbe)
                                } 
                                else {
                                        volcanoPlot(fit_tcx)
                                }
                        })
                        output$PCA <- renderPlotly({
                                # if(is.null(pcaEnv$PCA_cbe)){
                                #         print("pcaEnv is null")
                                # } else {
                                #         print("pcaEnv is not null")
                                # }
                                
                                if (region() == "CBE") {
                                        pcaPlot(pca_cbe, df_diag_cbe)
                                } else {
                                        pcaPlot(pca_tcx, df_diag_tcx)
                                }
                        })
                }
        )
}