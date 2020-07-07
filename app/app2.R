source("de.R")
source("pca.R")
source("pathway.R")

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
                radioButtons("region", "Brain Regions",
                             choices = c("Cerebellum" = "CBE",
                                         "Temporal cortex" = "TCX"),
                             selected = "CBE"
                )
        )
)

body <- dashboardBody(
        tags$link(href="styles.css", rel="stylesheet"),
        tabItems(
                # the first subtab content
                tabItem(
                        tabName = "rna-de",
                        deUI("rna-de")
                ),
                # the second subtab content
                tabItem(
                        tabName = "rna-pca",
                        pcaUI("rna-pca")
                ),
                # the third subtab content
                tabItem(
                        tabName = "rna-pathway",
                        pathwayUI("rna-pathway")
                )
        )
)

ui <- dashboardPage(
        dashboardHeader(title = "MayoRNAseq"),
        sidebar, body
)

server <- function(input, output, session) {
        deServer("rna-de", region = reactive(input$region))
        pcaServer("rna-pca", region = reactive(input$region))
        pathwayServer("rna-pathway", region = reactive(input$region))
}

shinyApp(ui, server)

