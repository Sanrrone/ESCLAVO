rm(list = ls())
library(shiny)
library(shinyDirectoryInput) #devtools::install_github('wleepang/shiny-directory-input')
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
#library(dashboardthemes) #install_github("nik01010/dashboardthemes")
library(base64enc)
library(filesstrings)
library(DT)
#library(jsonlite)
#source("Modules/themes_modules.R")
source("Classes/generic_pipeline_structure.R")
source("Modules/functions.R")
source("Modules/header_module.R")
source("Modules/rsidebar_ui.R")
source("Modules/analysisAvails_module.R")
source("Modules/callTabDef_module.R")
source("Modules/callTabContent_module.R")
source("Modules/analysisAvails_module.R")



#analysisAvail=c("16s18sits"=1,"rnaseq"=2,"metagenomic"=3,"gbs"=4)
#source from functions.R


dashboardPagePlus(
  #"blue", "blue-light", "black", "black-light", "purple", "purple-light", 
  #"green", "green-light", "red", "red-light", "yellow", "yellow-light"
  skin = "red-light",
  header = headerUIm("headermodule"),
  sidebar = dashboardSidebar(
    callTabDefUIm("tabdefmodule")
  ),
  rightsidebar = rightsidebarUIm("rsidebarmodule"),
  
  body = dashboardBody(
   #uiChangeThemeOutput(),
    tags$head(tags$style(HTML('
      .main-header .logo {
        font-family: "Georgia", Times, "Times New Roman", serif;
        font-weight: bold;
        font-size: 25px;
      }'
    ))),
   tags$head(tags$style(".modal-content {background: transparent;
                         box-shadow: inset 1px 2000px rgba(208, 208, 208, 0)}")),
   tags$head(
     tags$style(".content-wrapper { min-width: 800px !important; }")
   ),

   tabItems(
    callTabContent("tabcontmodule")
   )
  ),
  title = "ESCLAVO-PROJECT"
)

