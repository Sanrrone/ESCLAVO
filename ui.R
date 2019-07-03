library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(base64enc)
source("Modules/header_ui.R")
source("Modules/rsidebar_ui.R")
source("Modules/functions.R")

#analysisAvail=c("16s18sits"=1,"rnaseq"=2,"metagenomic"=3,"gbs"=4)
#source from functions.R

dashboardPagePlus(
  skin = "blue-light",
  header = headerUIm("headerUImodule",analysisAvail),
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("welcome", tabName = "Welcome", icon = icon("dashboard"),selected = T),
      uiOutput("tabsDef")
    )
  ),
  rightsidebar = uiOutput("rsidebarcontent"),
  
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
                  uiOutput("itemListUI")
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
              
      )
    )
  ),
  title = "ESCLAVO-PROJECT"
)

