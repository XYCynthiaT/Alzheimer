

glycotbl <- function(region, inOther = FALSE, class="N-glycosylation"){
    AD <- DE[[region]]$AD
    MCI <- DE[[region]]$MCI
    Other <- DE[[region]]$Other
    anova <- DE[[region]]$anova
    
    # DE of all glycosylation related genes
    AD_glyco <- AD %>%
        dplyr::filter(ensembl_gene_id %in% reactome$ensembl_gene_id) %>%
        right_join(distinct(reactome[,4:5]), by = "ensembl_gene_id") %>%
        dplyr::select(1:4, 9, 5:8)
    AD_glyco <- AD_glyco[complete.cases(AD_glyco),]
    # Exlude genes chaged in other dementia
    sigAD <- ifelse(AD$padj>0.05, 0, ifelse(AD$logFC<0, -1, 1))
    sigOther <- ifelse(Other$padj>0.05, 0, ifelse(Other$logFC<0, -1, 1))
    sigBoth <- sigAD == sigOther
    keep <- (sigAD != 0)&!sigBoth
    # Reactome genes
    glycoGenes <- AD$ensembl_gene_id %in% reactome$ensembl_gene_id
    AD_glyco_trend <- AD %>%
        dplyr::filter(keep&glycoGenes) %>%
        right_join(distinct(reactome[,4:5]), by = "ensembl_gene_id") %>%
        dplyr::filter(complete.cases(.)) %>%
        dplyr::select(1:4, 9, 5:8) 
    
    MCI_glyco <- MCI %>%
        right_join(AD_glyco[,c(1,2,5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
        dplyr::select(1:4, 9, 5:8)
    MCI_glyco_trend <- MCI_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_trend$ensembl_gene_id)
    
    Other_glyco <- Other %>%
        right_join(AD_glyco[,c(1,2,5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
        dplyr::select(1:4, 9, 5:8)
    Other_glyco_trend <- Other_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_trend$ensembl_gene_id)
    
    anova_glyco <- anova %>%
        right_join(AD_glyco[,c(1,2,5)], by = c("ensembl_gene_id", "hgnc_symbol")) %>%
        dplyr::select(1:4, 8, 5:7)
    anova_glyco_trend <- anova_glyco %>%
        dplyr::filter(ensembl_gene_id %in% AD_glyco_trend$ensembl_gene_id)
    if (inOther) {
        tbl <- left_join(AD_glyco, MCI_glyco[, c(1, 5:9)], by=c("ensembl_gene_id", "Class"), suffix = c(".AD", ".MCI")) %>%
            left_join(Other_glyco[, c(1, 5:9)], by=c("ensembl_gene_id", "Class")) %>%
            dplyr::rename(logFC.Other=logFC) %>%
            left_join(anova_glyco[, c(1, 5:8)], by=c("ensembl_gene_id", "Class"), suffix = c(".Other", ".ANOVA"))
    } else {
        tbl <- left_join(AD_glyco_trend, MCI_glyco_trend[, c(1, 5:9)], by=c("ensembl_gene_id", "Class"), suffix = c(".AD", ".MCI")) %>%
            left_join(Other_glyco_trend[, c(1, 5:9)], by=c("ensembl_gene_id", "Class")) %>%
            dplyr::rename(logFC.Other=logFC) %>%
            left_join(anova_glyco_trend[, c(1, 5:8)], by=c("ensembl_gene_id", "Class"), suffix = c(".Other", ".ANOVA"))
    }
    if (class!="all") {
        tbl <- filter(tbl, Class==class)
    } else {
        AD[, 5:ncol(AD)] <- signif(AD[, 5:ncol(AD)], 2)
        AD
    }
}
boxplotOutput <- function(tbl, gene_selected, region){
    rownames(DE[[1]]$AD) <- DE[[1]]$AD$ensembl_gene_id
    rownames(DE[[1]]$MCI) <- DE[[1]]$MCI$ensembl_gene_id
    rownames(DE[[1]]$Other) <- DE[[1]]$Other$ensembl_gene_id
    rownames(DE[[1]]$anova) <- DE[[1]]$anova$ensembl_gene_id
    ens <- tbl$ensembl_gene_id[gene_selected]
    df <- data.frame(
        count = logcpm$edata[ens,],
        diagnosis = logcpm$pdata$diagnosis
    ) 
    stat.test <- tibble(
        group1 = rep(1, 3),
        group2 = c(2, 3, 4),
        p.adj = c(DE[[region]]$MCI[ens,'padj'],DE[[region]]$AD[ens,'padj'], DE[[region]]$Other[ens,'padj']) %>%
            signif(digits = 3)
    ) 
    ggplot(df, aes(diagnosis, count, color = diagnosis)) +
        geom_boxplot()+
        geom_jitter(width = 0.1) +
        ylab("logcpm") +
        labs(title = paste0(tbl$hgnc_symbol[gene_selected], " in ", region, ", ANOVA Test: p.adj=", signif(DE[[region]]$anova[ens,'padj'], digits = 3))) +
        theme_bw() +
        stat_pvalue_manual(
            stat.test, 
            y.position = max(df$count), step.increase = 0.1,
            label = "p.adj"
        )
}
volcanoOutput_rosmap <- function(tbl){
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
rosmapUI <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            inputPanel(
                selectInput(ns("region"), 
                            "Brain Region", 
                            choices = c("Dorsolateral Prefrontal Cortex (DLPFC)"="dorsolateral")),
                radioButtons(ns("inother"), "Changed in other dementia too?", choices = c(TRUE, FALSE), inline = T),
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

rosmapServer <- function(id, region) {
    moduleServer(
        id,
        function(input, output, session){
            output$de <- renderDT({
                tbl <- glycotbl(input$region, input$inother, input$class)
                datatable(tbl, selection = list(mode="single", selected=1))%>%
                    formatSignif(columns = 6:ncol(tbl))
            })
            output$boxplot <- renderPlot({
                tbl <- glycotbl(input$region, input$inother, input$class)
                boxplotOutput(tbl, input$de_rows_selected, input$region)
            })
            output$vol <- renderPlotly({
                tbl <- glycotbl(input$region, input$inother, input$class)
                volcanoOutput_rosmap(tbl)
            })
        }
    )
}
# server = function(input, output, session) {
#         rosmapServer("a")
# }
# 
# shinyApp(ui = rosmapUI("a"), server = server)