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
                      img(src = "https://image.flaticon.com/icons/svg/204/204074.svg"), 
                      user_name = "16S", 
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

