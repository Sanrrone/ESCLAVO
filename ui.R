library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyBS)
library(base64enc)
source("Classes/generic_pipeline_structure.R")
source("Modules/functions.R")
source("Modules/header_ui.R")
source("Modules/rsidebar_ui.R")
source("Modules/analysisAvails_ui.R")
source("Modules/mainButtons_ui.R")
source("Modules/news_ui.R")


#analysisAvail=c("16s18sits"=1,"rnaseq"=2,"metagenomic"=3,"gbs"=4)
#source from functions.R

dashboardPagePlus(
  skin = "red-light",
  header = headerUIm("headerUImodule",analysisAvail),
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("welcome", tabName = "Welcome", icon = icon("dashboard"),selected = T),
      uiOutput("tabsDef")
    )
  ),
  rightsidebar = rightsidebarUIm("rsidebarmodule",analysisAvail),
  
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
                analysisAvailUIm("analysisavails",analysisAvail)
              ),
              mainButtonUIm("mainButtons"),
              newsUIm("news")
      )
    )
  ),
  title = "ESCLAVO-PROJECT"
)

