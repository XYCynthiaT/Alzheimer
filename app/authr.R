authrUI <- function(id) {
    ns <- NS(id)
    tagList(
        fluidPage(
            shinyjs::useShinyjs(),
            loginUI(ns("login")) # This is the inner UI
        )
    )
}

authrServer <- function(input, output, session) {
    # This is the inner server
    # call login module supplying data frame, user and password cols
    # and reactive trigger
    credentials <- callModule(shinyauthr::login,
                   id = "login",
                   data = user_base,
                   user_col = user,
                   pwd_col = password,
                   log_out = reactive(logout_init()))

    # call the logout module with reactive trigger to hide/show
    logout_init <- callModule(
        shinyauthr::logout, 
        "logout", 
        reactive(credentials()$user_auth)
    )
    # pulls out the user information returned from login module
    user_data <- reactive({credentials()$info})
    return(credentials)
}

# ui <- authrUI("authr")
# server <- function(input, output, session) {
#     callModule(authrServer, "authr")
# }
# shinyApp(ui, server)