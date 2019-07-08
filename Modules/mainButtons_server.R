source("Modules/startProject_server.R")

mainButtonsModule<-function(input, output, session) {
  ns <- session$ns
  
callModule(newProjectModule,"startProjectmodule")
#callModule(openProjectModule,"openProjectmodule")
#callModule(changeProjectDirModule,"changeProjectDirmodule")

}