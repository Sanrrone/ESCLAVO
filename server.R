source("Modules/rsidebar_server.R")
source("Modules/startProject_server.R")

function(input, output, session) {

  #callModule(module = serverChangeTheme, id = "moduleChangeTheme")
  callModule(headerModule,"headermodule")
  callModule(rsidebarModule,"rsidebarmodule")
  callModule(tabDefModule,"tabdefmodule")
  callModule(analysisAvailsModule,"analysisavails")
  callModule(newProjectModule,"startProjectmodule")
  #callModule(openProjectModule,"openProjectmodule")
  #callModule(changeProjectDirModule,"changeProjectDirmodule")
  callModule(tabContentModule,"tabcontmodule",session)

}
