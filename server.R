source("Modules/rsidebar_server.R")
source("Modules/callTabDef_server.R")
source("Modules/startProject_server.R")
source("Modules/analysisAvails_server.R")
source("Modules/callTabContent_server.R")

function(input, output) {
  
  callModule(rsidebarModule,"rsidebarmodule")
  callModule(tabDefModule,"tabdefmodule")
  callModule(analysisAvailsModule,"analysisavails")
  callModule(newProjectModule,"startProjectmodule")
  #callModule(openProjectModule,"openProjectmodule")
  #callModule(changeProjectDirModule,"changeProjectDirmodule")
  callModule(tabContentModule,"tabcontmodule")
  


}
