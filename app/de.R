# pkgs <- c("shiny", "dplyr", "ggplot2", "plotly", "shinydashboard", "HTSet", "DT")
# for (pkg in pkgs) {
#         suppressPackageStartupMessages(library(pkg, character.only = T))
# }

rna <- readRDS("data/RNA_normalized.RDS")
load("data/differentialExp.rda")

boxplotOutput <- function(gene_selected, region) {
        rest <- c("CBE", "TCX")[!(c("CBE", "TCX") %in% region)]
        if ((gene_selected %in% rownames(rna[[rest]]$edata))) {
                gene_cbe <- rna$CBE[gene_selected]
                gene_tcx <- rna$TCX[gene_selected]
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
                        geom_boxplot() +
                        geom_jitter(alpha = 0.3) +
                        labs(title = gene_selected) +
                        facet_wrap(~tissue) +
                        theme_bw()
        } else {
                gene <- rna[[region]][gene_selected]
                df <- data.frame(
                        expression = gene$edata[1,],
                        diagnosis = gene$pdata$Diagnosis
                )
                ggplot(df, aes(diagnosis, expression, color = diagnosis)) +
                        geom_boxplot() +
                        geom_jitter(alpha = 0.3) +
                        labs(title = paste0(gene_selected, ", ", region)
                        ) +
                        theme_bw()
        }
}

deUI <- function(id) {
        ns <- NS(id)
        tagList(
                column(
                        6,
                        box(
                                title = "Differential Expression",
                                width = NULL,
                                solidHeader = T,
                                status = "primary",
                                dataTableOutput(ns("de"))
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
                                plotlyOutput(ns("boxplot")),
                                textOutput(ns("txt"))
                        )
                )
        )
}

deServer <- function(id, region) {
        moduleServer(
                id,
                function(input, output, session) {
                        output$de <- renderDataTable({
                                if (region() == 'CBE') {
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
                        output$boxplot <- renderPlotly(
                                if (region() == 'CBE') {
                                        gene_selected <- output_lm_cbe$ensembl_id[input$de_rows_selected]
                                        boxplotOutput(gene_selected, region())
                                } else {
                                        gene_selected <- output_lm_tcx$ensembl_id[input$de_rows_selected]
                                        boxplotOutput(gene_selected, region())
                                }
                        )
                        
                        output$txt <- renderText({
                                if (region() == 'CBE') {
                                        gene_selected <- output_lm_cbe$ensembl_id[input$de_rows_selected]
                                        if (gene_selected %in% featureNames(rna$TCX) == F){
                                                return(
                                                        paste0(
                                                                "Notes: ",
                                                                gene_selected,
                                                                " is less than 10 count per million in the temporal cortex."
                                                        )
                                                )
                                                
                                        }
                                } else {
                                        gene_selected <- output_lm_tcx$ensembl_id[input$de_rows_selected]
                                        if (gene_selected %in% featureNames(rna$CBE) == F) {
                                                return(
                                                        paste0(
                                                                "Notes: ",
                                                                gene_selected,
                                                                " is less than 10 count per million in the cerebellum."
                                                        )
                                                )
                                        }
                                        
                                } 
                        })
                }
        )
}
