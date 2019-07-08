source("Modules/startProject_ui.R")
source("Modules/openProject_ui.R")
source("Modules/changeProjectDir_ui.R")

mainButtonsUIm<-function(id){
  ns<-NS(id)
  
  fluidRow(
    newProjectUIm(ns("startProjectmodule")),
    openProjectUIm(ns("openProjectmodule")),
    changeProjectDirUIm(ns("changeProjectDirmodule"))
  )
}
