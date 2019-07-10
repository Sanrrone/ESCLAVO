source("Modules/startProject_ui.R")
source("Modules/openProject_ui.R")
source("Modules/changeProjectDir_ui.R")
source("Modules/news_ui.R")
source("Modules/analysisAvails_ui.R")

callTabContent<-function(id){
  ns<-NS(id)
  
  
  uiOutput(ns("tabContentUI"))
  
}