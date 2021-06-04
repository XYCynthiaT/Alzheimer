source("global.R")
source("mayo.R")
source("msbb.R")
source("rosmap.R")
source("regions.R")

dashboardUI <- function(id) {
        ns <- NS(id)
        header <- dashboardHeader(
                title = "AD RNAseq",
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
                        menuItem("Mayo (TC, CBE)", tabName = "mayo"),
                        # second menu
                        menuItem("MSBB (FP, STG, IFG, PHG)", tabName = "msbb"),
                        # third menu
                        menuItem("ROSMAP (DLPFC)", tabName = "rosmap"),
                        menuItem("Summary across regions", tabName = "regions")
                )
        )
        
        body <- dashboardBody(
                shinyjs::useShinyjs(),
                tags$link(href="styles.css", rel="stylesheet"),
                tabItems(
                        # the first subtab content
                        tabItem(
                                tabName = "mayo",
                                mayoUI(ns("mayo"))
                        ),
                        # the second subtab content
                        tabItem(
                                tabName = "msbb",
                                msbbUI(ns("msbb"))
                        ),
                        # the third subtab content
                        tabItem(
                                tabName = "rosmap",
                                rosmapUI(ns("rosmap"))
                        ),
                        tabItem(
                                tabName = "regions",
                                regsUI(ns("regions"))
                        )
                )
        )
        dashboardPage(header, sidebar, body)
}
        

dashboardServer <- function(id) {
        moduleServer(
                id,
                function(input, output, session) {
                        mayoServer("mayo")
                        msbbServer("msbb")
                        rosmapServer("rosmap")
                        regsServer("regions")
                }
        )
}
server = function(input, output, session) {
        dashboardServer("a")
}

shinyApp(ui = dashboardUI("a"), server = server)