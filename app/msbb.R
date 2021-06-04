
glycotbl_msbb <- function(region, class="N-glycosylation"){
    AD <- DE_msbb[[region]]$AD

    # DE of all glycosylation related genes
    AD_glyco <- AD %>%
        dplyr::filter(ensembl_gene_id %in% reactome$ensembl_gene_id) %>%
        right_join(distinct(reactome[,4:5]), by = "ensembl_gene_id") %>%
        dplyr::select(1:4, 9, 5:8)
    AD_glyco <- AD_glyco[complete.cases(AD_glyco),]
    if (class!="all") {
        tbl <- AD_glyco %>%
            filter(Class==class)
    } else {
        AD[, 5:ncol(AD)] <- signif(AD[, 5:ncol(AD)], 2)
        AD
    }
}
boxplotOutput_msbb <- function(tbl, gene_selected, region){
    rownames(DE_msbb[[region]]$AD) <- DE_msbb[[region]]$AD$ensembl_gene_id
    ens <- tbl$ensembl_gene_id[gene_selected]
    df <- data.frame(
        count = logcpm_msbb[[region]]$edata[ens,],
        diagnosis = logcpm_msbb[[region]]$pdata$diagnosis
    )  %>%
        filter(diagnosis != "Transition")
    stat.test <- tibble(
        group1 = 1,
        group2 = c(2),
        p.adj = c(DE_msbb[[region]][["AD"]][ens,'padj']) %>%
            signif(digits = 3)
    )  
    ggplot(df, aes(diagnosis, count, color = diagnosis)) +
        geom_boxplot()+
        geom_jitter(width = 0.1) +
        ylab("logcpm") +
        labs(title = paste0(tbl$hgnc_symbol[gene_selected], " in ", region)) +
        theme_bw() +
        stat_pvalue_manual(
            stat.test, 
            y.position = max(df$count), step.increase = 0.1,
            label = "p.adj"
        )
}
volcanoOutput_msbb <- function(tbl){
    tbl %>%
        mutate(sign=ifelse(pval>0.05, 0, ifelse(logFC>0, 1, -1))) %>%
        mutate(sign=factor(sign, levels = c(0, -1, 1), labels = c("Not sig", "Decrease in AD", "Increase in AD"))) %>%
        ggplot(aes(logFC, -log(pval)))+
        geom_point(aes(color=sign, symbol=hgnc_symbol, ensembl=ensembl_gene_id))+
        geom_hline(yintercept = -log(0.05))+
        scale_color_manual(values = c("grey", "steelblue", "red"))+
        theme_bw()
    
}

msbbUI <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            inputPanel(
                selectInput(ns("region"), 
                             "Brain Region", 
                             choices = c("Frontal Pole(FP)"="BM10", "Inferior Frontal Gyrus (IFG)"="BM44", 
                                         "Superior Temporal Gyrus (STG)"="BM22", "Parahippocampal Gyrus (PHG)"="BM36"), 
                             selected = "BM10"),
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

msbbServer <- function(id, region) {
    moduleServer(
        id,
        function(input, output, session){
            output$de <- renderDT({
                tbl <- glycotbl_msbb(input$region, input$class)
                datatable(tbl, selection = list(mode="single", selected=1))%>%
                    formatSignif(columns = 6:ncol(tbl))
            })
            output$boxplot <- renderPlot({
                tbl <- glycotbl_msbb(input$region, input$class)
                boxplotOutput_msbb(tbl, input$de_rows_selected, input$region)
            })
            output$vol <- renderPlotly({
                tbl <- glycotbl_msbb(input$region, input$class)
                volcanoOutput_msbb(tbl)
            })
        }
    )
}
# server = function(input, output, session) {
#         msbbServer("a")
# }
# 
# shinyApp(ui = msbbUI("a"), server = server)