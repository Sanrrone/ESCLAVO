callTabDefUIm<-function(id, projectName){
  ns <- NS(id)
  sidebarMenuOutput(ns("tabsDefUI"))
  
}


tabDefModule<-function(input, output, session){
  ns <- session$ns

  output$tabsDefUI<- renderMenu({

    if(projectName()!=""){
      sidebarMenu(id = "mainMenu",
              menuItem("welcome", tabName = "Welcome", icon = icon("dashboard")),
              menuItem(paste0("Project: ",projectName()), tabName = "projectStatus",
                            icon = icon("chalkboard-teacher"),selected = T),
              lapply(getAnalysisSteps(pipelines[[analysisType()]]),function(x){
                menuItem(x$stepName, tabName = paste0("step_",x$stepID), icon = x$icon)
              }),
              menuItem("Static/Dynamic Report", tabName = "projectReport", icon = icon("file-invoice"))
        )
      
    }else{
      sidebarMenu(menuItem("welcome", tabName = "Welcome", icon = icon("dashboard")))
    }
  })
  
}
