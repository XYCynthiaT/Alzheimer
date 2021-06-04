box_regs <- function(symbol="LRP1", ADonly=FALSE){
    ens <- featureNames(logcpm)[logcpm$fdata$hgnc_symbol==symbol]
    dlpfc <- data.frame(
        Region = "DLPFC",
        Diagnosis = logcpm$pdata$diagnosis,
        Logcpm = logcpm$edata[ens,]
    )
    tcx <- data.frame(
        Region = "TCX",
        Diagnosis = logcpm_mayo$TCX$pdata$Diagnosis,
        Logcpm = logcpm_mayo$TCX$edata[ens,]
    ) %>%
        mutate(Diagnosis = sub("Con", "NCI", Diagnosis))
    df <- rbind(dlpfc, tcx)
    for (reg in names(logcpm_msbb)) {
        bm <- data.frame(
            Region = reg,
            Diagnosis = logcpm_msbb[[reg]]$pdata$diagnosis,
            Logcpm = logcpm_msbb[[reg]]$edata[ens,]
        )
        df <- rbind(df,bm)
    }
    cbe <- data.frame(
        Region = "CBE",
        Diagnosis = logcpm_mayo$CBE$pdata$Diagnosis,
        Logcpm = logcpm_mayo$CBE$edata[ens,]
    ) %>%
        mutate(Diagnosis = sub("Con", "NCI", Diagnosis))
    df <- rbind(df, cbe)
    if (ADonly) {
        df <- df %>%
            filter(Diagnosis %in% c('NCI', "AD")) %>%
            mutate(Diagnosis = factor(Diagnosis, levels = c("NCI", "AD")))
    } else {
        df <- df %>%
            filter(Diagnosis != "Transition") %>%
            mutate(Diagnosis = factor(Diagnosis, levels = c("NCI", "MCI", "PA", "AD", "PSP", "Other")))
        
    }
    # rownames(DE_mayo[[region]]$AD) <- DE_mayo[[region]]$AD$ensembl_gene_id
    # rownames(DE_mayo[[region]]$PA) <- DE_mayo[[region]]$AD$ensembl_gene_id
    # rownames(DE_mayo[[region]]$PSP) <- DE_mayo[[region]]$AD$ensembl_gene_id
    # rownames(DE_mayo[[region]]$anova) <- DE_mayo[[region]]$AD$ensembl_gene_id
    # stat.test <- tibble(
    #     group1 = rep(6.6, 3),
    #     group2 = c(6.8, 7.2, 7.4),
    #     p.adj = c(DE_mayo[[region]]$PA[ens,'padj'], DE_mayo[[region]]$PSP[ens,'padj'], DE_mayo[[region]]$AD[ens,'padj']) %>%
    #         signif(digits = 3)
    # )
    ggplot(df, aes(Region, Logcpm, color=Diagnosis))+
        geom_boxplot(aes(color=Diagnosis))+
        geom_point(position = position_jitterdodge(jitter.width = 0.1))+
        theme_bw() 
        # stat_pvalue_manual(
        #     stat.test, 
        #     y.position = max(df$Logcpm), step.increase = 0.1,
        #     label = "p.adj"
        # )
}
pval <- function(symbol="LRP1"){
    ens <- featureNames(logcpm)[logcpm$fdata$hgnc_symbol==symbol]
    # mayo
    mayo1 <- lapply(c("AD", "PA", "PSP"), function(x){
        ind <- match(ens, DE_mayo[["TCX"]][[x]]$ensembl_gene_id)
        DE_mayo[["TCX"]][[x]][ind,]%>%
            mutate(Diagnosis=x, Region="TCX")
    }) 
    mayo2 <- lapply(c("AD", "PA", "PSP"), function(x){
        ind <- match(ens, DE_mayo[["CBE"]][[x]]$ensembl_gene_id)
        DE_mayo[["CBE"]][[x]][ind,]%>%
            mutate(Diagnosis=x, Region="CBE")
    }) 
    # ROSMAP
    rosmap1 <- lapply(c("AD", "MCI", "Other"), function(x){
        ind <- match(ens, DE$dorsolateral[[x]]$ensembl_gene_id)
        DE$dorsolateral[[x]][ind,]%>%
            mutate(Diagnosis=x, Region="DLPFC")
    }) 
    # MSBB
    msbb1 <- lapply(names(DE_msbb), function(x){
        ind <- match(ens, DE_msbb[[x]][["AD"]]$ensembl_gene_id)
        DE_msbb[[x]][["AD"]][ind,]%>%
            mutate(Diagnosis="AD", Region=x)
    })
    list(mayo1, mayo2, rosmap1, msbb1)%>%
        Reduce(append, .) %>%
        Reduce(rbind, .) %>%
        select(Region, Diagnosis, everything())
}

regsUI <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            inputPanel(
                h4("Input a hgnc_symbol of a gene. Make sure all characters are in capital."),
                textInput(ns("symbol"), "hgnc_symbol"),
                actionButton(ns("ok"), "OK")
            )
        ),
        fluidRow(
            column(12, textOutput(ns("check")))
        ),
        fluidRow(
            box(
                title = "Gene Expression across the brain (NCI vs AD)",
                width = 12,
                solidHeader = T,
                status = "primary",
                plotOutput(ns("box_AD"))
            )
        ),
        fluidRow(
            box(
                title = "Gene Expression across the brain (all Dx)",
                width = 12,
                solidHeader = T,
                status = "primary",
                plotOutput(ns("box_allDx"))
            )
        ),
        fluidRow(
            box(
                title = "Statistical Results",
                width = 12,
                solidHeader = T,
                status = "primary",
                DTOutput(ns("pval"))
            )
        )
    )
}

regsServer <- function(id, region) {
    moduleServer(
        id,
        function(input, output, session){
            check_text <- eventReactive(input$ok, {
                if(input$symbol %in% logcpm$fdata$hgnc_symbol){
                    paste0(input$symbol, " is found in RNA-seq data.")
                } else {
                    paste0(input$symbol, " is not found in RNA-seq data, please check the hgnc_symbol.")
                }
            })
            output$check <- renderText({
                check_text()
            })
            allDx <- eventReactive(input$ok, {
                if(input$symbol %in% logcpm$fdata$hgnc_symbol){
                    box_regs(input$symbol, ADonly=FALSE)
                }
            })
            AD_NCI <- eventReactive(input$ok, {
                if(input$symbol %in% logcpm$fdata$hgnc_symbol){
                    box_regs(input$symbol, ADonly=TRUE)
                }
            })
            pval_tbl <- eventReactive(input$ok, {
                if(input$symbol %in% logcpm$fdata$hgnc_symbol){
                    pval(symbol = input$symbol) %>%
                        datatable(rownames = FALSE) %>%
                        formatSignif(7:10)
                }
            })
            output$box_allDx <- renderPlot({
                allDx()
            })
            output$box_AD <- renderPlot({
                AD_NCI()
            })  
            output$pval <- renderDT({
                pval_tbl()
            }) 
        }
    )
}
# server = function(input, output, session) {
#         regsServer("a")
# }
# 
# shinyApp(ui = regsUI("a"), server = server)