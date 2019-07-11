library(shiny)
library(shinyDirectoryInput) #devtools::install_github('wleepang/shiny-directory-input')
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(dashboardthemes)
library(base64enc)
library(filesstrings)
#library(jsonlite)
source("Modules/themes_modules.R")
source("Classes/generic_pipeline_structure.R")
source("Modules/functions.R")
source("Modules/header_ui.R")
source("Modules/rsidebar_ui.R")
source("Modules/callTabDef_ui.R")
source("Modules/callTabContent_ui.R")


#analysisAvail=c("16s18sits"=1,"rnaseq"=2,"metagenomic"=3,"gbs"=4)
#source from functions.R


dashboardPagePlus(
  #skin = "black",
  header = headerUIm("headerUImodule"),
  sidebar = dashboardSidebar(
    callTabDefUIm("tabdefmodule")
  ),
  rightsidebar = rightsidebarUIm("rsidebarmodule"),
  
  body = dashboardBody(
   uiChangeThemeOutput(),
    tags$head(tags$style(HTML('
      .main-header .logo {
        font-family: "Georgia", Times, "Times New Roman", serif;
        font-weight: bold;
        font-size: 25px;
      }'
    ))),
    callTabContent("tabcontmodule")
    
  ),
  title = "ESCLAVO-PROJECT"
)

