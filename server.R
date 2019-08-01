source("Modules/rsidebar_server.R")
source("Modules/startProject_module.R")
source("Modules/openProject_module.R")

function(input, output, session) {
  
  #callModule(module = serverChangeTheme, id = "moduleChangeTheme")
  callModule(headerModule,"headermodule")
  callModule(rsidebarModule,"rsidebarmodule")
  callModule(tabDefModule,"tabdefmodule")
  callModule(analysisAvailsModule,"analysisavails")
  callModule(newProjectModule,"startProjectmodule")
  callModule(openProjectModule,"openProjectmodule")
  #callModule(changeProjectDirModule,"changeProjectDirmodule")
  callModule(tabContentModule,"tabcontmodule",session)
  
  #call Step pipelines modules
  observe({
    if(nrow(projectConf())!=0){
      for(x in getAnalysisSteps(pipelines[[analysisType()]])[1]){
        mVector<-x$tabcontentSrv
        callModule(module = mVector$server,id = mVector$id,x$stepID)
      }
    }
  })

  
}
