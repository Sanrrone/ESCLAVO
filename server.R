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
    if(projectName()!=""){
      
      callmodules<-reactivePoll(intervalMillis = 10000, session = session,
                                       checkFunc = function(){digest(paste0(projectFolder(),"/",projectName(),"_eConf.tsv"),
                                                                                       algo="md5",file=TRUE)},
                                       valueFunc = function(){
                                           lapply(getAnalysisSteps(pipelines[[analysisType()]]), function(x){
                                             mVector<-x$tabcontentSrv
                                             callModule(module = mVector$server,id = mVector$id, x$stepID)
                                           })
                             })
      callmodules()
    }
  })

  
}
