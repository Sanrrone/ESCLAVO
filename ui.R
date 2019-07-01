library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
source("Modules/header_ui.R")

dashboardPagePlus(
  header = headerUIm("headerUImodule"),
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("welcome", tabName = "Welcome", icon = icon("dashboard"),selected = T),
      uiOutput("tabsDef")
    )
  ),
  body = dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "Welcome",
              fluidRow(
                box(
                  title = "Available analysis",
                  status = "success",
                  width = 12,
                  userList(
                    userListItem(
                      src = "https://raw.githubusercontent.com/Sanrrone/ESCLAVO/master/images/16srrna.png", 
                      user_name = "16S", 
                      description = "28.04.2018"
                    ),
                    userListItem(
                      src = "https://www.rstudio.com/wp-content/uploads/2014/04/knitr.png", 
                      user_name = "knitr", 
                      description = "28.04.2018"
                    )
                  )
                )
              ),
              fluidRow(
                box(
                  title = "Create new project",
                  status = "success",
                  width = 6
                ),
                box(
                  title = "Open project",
                  status = "success",
                  width = 6
                )
              )
              
      ),
      uiOutput("tabs")
    )
  ),
  title = "ESCLAVO-PROJECT"
)

