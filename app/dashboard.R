dashboardUI <- function(id) {
        ns <- NS(id)
        header <- dashboardHeader(
                title = "MayoRNAseq",
                tags$li(
                        class = "nav-item",
                        tags$a(
                                class = "btn btn-danger action-button",
                                id = "logout-button",
                                type = "button",
                                icon("sign-out-alt"), "Log out"
                        )
                )
        )
        
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
                        radioButtons(ns("region"), "Brain Regions",
                                     choices = c("Cerebellum" = "CBE",
                                                 "Temporal cortex" = "TCX"),
                                     selected = "CBE"
                        )
                )
        )
        
        body <- dashboardBody(
                shinyjs::useShinyjs(),
                tags$link(href="styles.css", rel="stylesheet"),
                tabItems(
                        # the first subtab content
                        tabItem(
                                tabName = "rna-de",
                                deUI(ns("rna-de"))
                        ),
                        # the second subtab content
                        tabItem(
                                tabName = "rna-pca",
                                pcaUI(ns("rna-pca"))
                        ),
                        # the third subtab content
                        tabItem(
                                tabName = "rna-pathway",
                                pathwayUI(ns("rna-pathway"))
                        )
                )
        )
        dashboardPage(header, sidebar, body)
}
        

dashboardServer <- function(id) {
        moduleServer(
                id,
                function(input, output, session) {
                        deServer("rna-de", region = reactive(input$region))
                        pcaServer("rna-pca", region = reactive(input$region))
                        pathwayServer("rna-pathway", region = reactive(input$region))
                }
        )
}
# server = function(input, output, session) {
#         dashboardServer("a")
# }
# 
# shinyApp(ui = dashboardUI("a"), server = server)