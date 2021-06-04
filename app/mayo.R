# pkgs <- c("shiny", "dplyr", "ggplot2", "plotly", "shinydashboard", "HTSet", "DT")
# for (pkg in pkgs) {
#         suppressPackageStartupMessages(library(pkg, character.only = T))
# }
glycotbl_mayo <- function(region, inPSP = FALSE, class="N-glycosylation"){
        AD <- DE_mayo[[region]]$AD
        PA <- DE_mayo[[region]]$PA
        PSP <- DE_mayo[[region]]$PSP
        anova <- DE_mayo[[region]]$anova
        
        # DE of all-glycosylation related genes
        AD_glyco <- AD %>%
                dplyr::filter(ensembl_gene_id %in% reactome$ensembl_gene_id) %>%
                right_join(distinct(reactome[,4:5]), by = "ensembl_gene_id") %>%
                dplyr::select(1:4, 9, 5:8)
        AD_glyco <- AD_glyco[complete.cases(AD_glyco),]
        # Exlude genes chaged in PSP
        sigAD <- ifelse(AD$padj>0.05, 0, ifelse(AD$logFC<0, -1, 1))
        sigPSP <- ifelse(PSP$padj>0.05, 0, ifelse(PSP$logFC<0, -1, 1))
        sigBoth <- sigAD == sigPSP
        keep <- (sigAD != 0)&!sigBoth
        # Reactome genes
        glycoGenes <- AD$ensembl_gene_id %in% reactome$ensembl_gene_id
        AD_glyco_trend <- AD %>%
                dplyr::filter(keep&glycoGenes) %>%
                right_join(distinct(reactome[,4:5]), by = "ensembl_gene_id") %>%
                dplyr::filter(complete.cases(.)) %>%
                dplyr::select(1:4, 9, 5:8) 
        
        PA_glyco <- PA %>%
                right_join(AD_glyco[,c(1,2,5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
                dplyr::select(1:4, 9, 5:8)
        PA_glyco_trend <- PA_glyco %>%
                dplyr::filter(ensembl_gene_id %in% AD_glyco_trend$ensembl_gene_id)
        
        PSP_glyco <- PSP %>%
                right_join(AD_glyco[,c(1,2,5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
                dplyr::select(1:4, 9, 5:8)
        PSP_glyco_trend <- PSP_glyco %>%
                dplyr::filter(ensembl_gene_id %in% AD_glyco_trend$ensembl_gene_id)
        
        anova_glyco <- anova %>%
                right_join(AD_glyco[,c(1,2,5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
                dplyr::select(1:4, 8, 5:7)
        anova_glyco_trend <- anova_glyco %>%
                dplyr::filter(ensembl_gene_id %in% AD_glyco_trend$ensembl_gene_id)
        
        if (inPSP) {
                tbl <- left_join(AD_glyco, PA_glyco[, c(1, 5:9)], by=c("ensembl_gene_id", "Class"), suffix = c(".AD", ".PA")) %>%
                        left_join(PSP_glyco[, c(1, 5:9)], by=c("ensembl_gene_id", "Class")) %>%
                        dplyr::rename(logFC.PSP=logFC) %>%
                        left_join(anova_glyco[, c(1, 5:8)], by=c("ensembl_gene_id", "Class"), suffix = c(".PSP", ".ANOVA"))
        } else {
                tbl <- left_join(AD_glyco_trend, PA_glyco_trend[, c(1, 5:9)], by=c("ensembl_gene_id", "Class"), suffix = c(".AD", ".PA")) %>%
                        left_join(PSP_glyco_trend[, c(1, 5:9)], by=c("ensembl_gene_id", "Class")) %>%
                        dplyr::rename(logFC.PSP=logFC) %>%
                        left_join(anova_glyco_trend[, c(1, 5:8)], by=c("ensembl_gene_id", "Class"), suffix = c(".PSP", ".ANOVA"))
        }
        if (class!="all") {
                tbl <- filter(tbl, Class==class)
        } else {
                AD[, 5:ncol(AD)] <- signif(AD[, 5:ncol(AD)], 2)
                AD
        }
}
# boxplotOutput <- function(gene_selected, region) {
#         rest <- c("CBE", "TCX")[!(c("CBE", "TCX") %in% region)]
#         if ((gene_selected %in% rownames(rna[[rest]]$edata))) {
#                 gene_cbe <- rna$CBE[gene_selected]
#                 gene_tcx <- rna$TCX[gene_selected]
#                 df_cbe <- data.frame(
#                         expression = gene_cbe$edata[1,],
#                         diagnosis = gene_cbe$pdata$Diagnosis,
#                         tissue = gene_cbe$pdata$Tissue
#                 )
#                 df_tcx <- data.frame(
#                         expression = gene_tcx$edata[1,],
#                         diagnosis = gene_tcx$pdata$Diagnosis, 
#                         tissue = gene_tcx$pdata$Tissue
#                 )
#                 df <- rbind(df_cbe, df_tcx)
#                 ggplot(df, aes(diagnosis, expression, color = diagnosis)) +
#                         geom_boxplot() +
#                         geom_jitter(alpha = 0.3) +
#                         labs(title = gene_selected) +
#                         facet_wrap(~tissue) +
#                         theme_bw()
#         } else {
#                 gene <- rna[[region]][gene_selected]
#                 df <- data.frame(
#                         expression = gene$edata[1,],
#                         diagnosis = gene$pdata$Diagnosis
#                 )
#                 ggplot(df, aes(diagnosis, expression, color = diagnosis)) +
#                         geom_boxplot() +
#                         geom_jitter(alpha = 0.3) +
#                         labs(title = paste0(gene_selected, ", ", region)
#                         ) +
#                         theme_bw()
#         }
# }
boxplotOutput_mayo <- function(tbl, gene_selected, region){
        rownames(DE_mayo[[region]]$AD) <- DE_mayo[[region]]$AD$ensembl_gene_id
        rownames(DE_mayo[[region]]$PA) <- DE_mayo[[region]]$AD$ensembl_gene_id
        rownames(DE_mayo[[region]]$PSP) <- DE_mayo[[region]]$AD$ensembl_gene_id
        rownames(DE_mayo[[region]]$anova) <- DE_mayo[[region]]$AD$ensembl_gene_id
        ens <- tbl$ensembl_gene_id[gene_selected]
        df <- data.frame(
                count = logcpm_mayo[[region]]$edata[ens,],
                diagnosis = logcpm_mayo[[region]]$pdata$Diagnosis
        ) 
        stat.test <- tibble(
                group1 = rep(1, 3),
                group2 = c(2, 3, 4),
                p.adj = c(DE_mayo[[region]]$PA[ens,'padj'], DE_mayo[[region]]$PSP[ens,'padj'], DE_mayo[[region]]$AD[ens,'padj']) %>%
                        signif(digits = 3)
        ) 
        ggplot(df, aes(diagnosis, count, color = diagnosis)) +
                geom_boxplot()+
                geom_jitter(width = 0.1) +
                ylab("logcpm") +
                labs(title = paste0(tbl$hgnc_symbol[gene_selected], " in ", region, ", ANOVA Test: p.adj=", signif(DE_mayo[[region]]$anova[ens,'padj'], digits = 3))) +
                theme_bw() +
                stat_pvalue_manual(
                        stat.test, 
                        y.position = max(df$count), step.increase = 0.1,
                        label = "p.adj"
                )
}
volcanoOutput_mayo <- function(tbl){
        if ("logFC" %in% colnames(tbl)) {
                tbl %>%
                        mutate(sign=ifelse(pval>0.05, 0, ifelse(logFC>0, 1, -1))) %>%
                        mutate(sign=factor(sign, levels = c(0, -1, 1), labels = c("Not sig", "Decrease in AD", "Increase in AD"))) %>%
                        ggplot(aes(logFC, -log(pval)))+
                        geom_point(aes(color=sign, symbol=hgnc_symbol, ensembl=ensembl_gene_id), alpha=0.3)+
                        geom_hline(yintercept = -log(0.05))+
                        scale_color_manual(values = c("grey", "steelblue", "red"))+
                        theme_bw()
        } else {
                tbl %>%
                        mutate(sign=ifelse(pval.AD>0.05, 0, ifelse(logFC.AD>0, 1, -1))) %>%
                        mutate(sign=factor(sign, levels = c(0, -1, 1), labels = c("Not sig", "Decrease in AD", "Increase in AD"))) %>%
                ggplot(aes(logFC.AD, -log(pval.AD)))+
                        geom_point(aes(color=sign, symbol=hgnc_symbol, ensembl=ensembl_gene_id))+
                        geom_hline(yintercept = -log(0.05))+
                        scale_color_manual(values = c("grey", "steelblue", "red"))+
                        theme_bw()
        }
}

mayoUI <- function(id) {
        ns <- NS(id)
        tagList(
                fluidRow(
                        inputPanel(
                                selectInput(ns("region"), "Brain Region", choices = c("Temporal Cortex"="TCX", "Cerebellum"="CBE"), selected = "TCX"),
                                radioButtons(ns("inPSP"), "Changed in PSP too?", choices = c(TRUE, FALSE), inline = T),
                                selectInput(ns("class"), "Glycosylation Class", choices = c("all", "N-glycosylation", "O-glycosylation", "Glycosphingolipid"), selected = "N-glycosylation")
                        )
                ),
                fluidRow(
                        box(
                                title = "Glycosylation Genes",
                                width = 6,
                                solidHeader = T,
                                status = "primary",
                                dataTableOutput(ns("de"))
                        ),
                        box(
                                title = "Boxplot",
                                width = 6,
                                solidHeader = T,
                                status = "primary",
                                plotOutput(ns("boxplot"))
                        )
                ),
                fluidRow(
                        column(6, offset = 6,
                               box(
                                       title = "Volcano Plot",
                                       width = NULL,
                                       solidHeader = T,
                                       status = "primary",
                                       plotlyOutput(ns("vol"))
                               ))
                )
        )
}

mayoServer <- function(id, region) {
        moduleServer(
                id,
                # function(input, output, session) {
                #         output$de <- renderDataTable({
                #                 if (region() == 'CBE') {
                #                         datatable(output_lm_cbe, 
                #                                   rownames = FALSE,
                #                                   selection = list(mode = 'single', selected = 1)
                #                         ) %>%
                #                                 formatSignif(columns = 2:7, digits = 2)
                #                 } else {
                #                         datatable(output_lm_tcx, 
                #                                   rownames = FALSE,
                #                                   selection = list(mode = 'single', selected = 1)
                #                         ) %>%
                #                                 formatSignif(columns = 2:7, digits = 2)
                #                 }
                #         })
                #         output$boxplot <- renderPlotly(
                #                 if (region() == 'CBE') {
                #                         gene_selected <- output_lm_cbe$ensembl_id[input$de_rows_selected]
                #                         boxplotOutput(gene_selected, region())
                #                 } else {
                #                         gene_selected <- output_lm_tcx$ensembl_id[input$de_rows_selected]
                #                         boxplotOutput(gene_selected, region())
                #                 }
                #         )
                #         
                #         output$txt <- renderText({
                #                 if (region() == 'CBE') {
                #                         gene_selected <- output_lm_cbe$ensembl_id[input$de_rows_selected]
                #                         if (gene_selected %in% featureNames(rna$TCX) == F){
                #                                 return(
                #                                         paste0(
                #                                                 "Notes: ",
                #                                                 gene_selected,
                #                                                 " is less than 10 count per million in the temporal cortex."
                #                                         )
                #                                 )
                #                                 
                #                         }
                #                 } else {
                #                         gene_selected <- output_lm_tcx$ensembl_id[input$de_rows_selected]
                #                         if (gene_selected %in% featureNames(rna$CBE) == F) {
                #                                 return(
                #                                         paste0(
                #                                                 "Notes: ",
                #                                                 gene_selected,
                #                                                 " is less than 10 count per million in the cerebellum."
                #                                         )
                #                                 )
                #                         }
                #                         
                #                 } 
                #         })
                # }
                function(input, output, session){
                        output$de <- renderDT({
                                tbl <- glycotbl_mayo(input$region, as.logical(input$inPSP), input$class)
                                datatable(tbl, selection = list(mode="single", selected=1))%>%
                                        formatSignif(columns = 6:ncol(tbl))
                        })
                        output$boxplot <- renderPlot({
                                tbl <- glycotbl_mayo(input$region, as.logical(input$inPSP), input$class)
                                boxplotOutput_mayo(tbl, input$de_rows_selected, input$region)
                        })
                        output$vol <- renderPlotly({
                                tbl <- glycotbl_mayo(input$region, as.logical(input$inPSP), input$class)
                                volcanoOutput_mayo(tbl)
                        })
                }
        )
}
# server = function(input, output, session) {
#         mayoServer("a")
# }
# 
# shinyApp(ui = mayoUI("a"), server = server)