library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(base64enc)
source("Modules/header_ui.R")

dashboardPagePlus(
  skin = "blue-light",
  header = headerUIm("headerUImodule"),
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("welcome", tabName = "Welcome", icon = icon("dashboard"),selected = T),
      uiOutput("tabsDef")
    )
  ),
  body = dashboardBody(
    tags$head(tags$style(HTML('
      .main-header .logo {
        font-family: "Georgia", Times, "Times New Roman", serif;
        font-weight: bold;
        font-size: 25px;
      }'
    ))),
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
                      src = dataURI(file="images/16srrna.png", mime="image/png"), 
                      user_name = "16S, 18S/ITS", 
                      description = "Updated: August 2019"
                    ),
                    userListItem(
                        src = dataURI(file="images/rnaseq.jpg", mime="image/png"), 
                        user_name = "RNA-seq", 
                        description = "Updated: August 2019"
                    ),
                    userListItem(
                      src = dataURI(file="images/metagenomic.jpg", mime="image/png"), 
                      user_name = "Metagenomic", 
                      description = "Updated: August 2019"
                    ),
                    userListItem(
                      src = dataURI(file="images/gbs.jpg", mime="image/png"), 
                      user_name = "Genotype by sequencing", 
                      description = "Updated: August 2019"
                    )
                  )
                )
              ),
              fluidRow(
                box(
                  title = "Create new project",
                  status = "success",
                  width = 4,
                  actionBttn(
                    inputId = "newprojectbtn",
                    label = "Start",
                    style = "material-flat",
                    color = "success",
                    icon = icon("rocket")
                  )
                ),
                box(
                  title = "Search existing project",
                  status = "success",
                  width = 4,
                  actionBttn(
                    inputId = "openprojectbtn",
                    label = "Open",
                    style = "material-flat",
                    color = "success",
                    icon = icon("folder")
                  )
                  #shinyDirButton("openprojectbtn", "Chose directory", "Upload")
                ),
                box(
                  title = "Change project directory",
                  status = "success",
                  width = 4,
                  actionBttn(
                    inputId = "changeprojectdirbtn",
                    label = "Choose",
                    style = "material-flat",
                    color = "success",
                    icon = icon("gears")
                  )
                )
              )
              
      ),
      uiOutput("tabs")
    )
  ),
  title = "ESCLAVO-PROJECT"
)

