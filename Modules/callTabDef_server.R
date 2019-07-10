tabDefModule<-function(input, output, session){
  ns <- session$ns

  output$tabsDefUI<- renderMenu({

    if(projectName()!=""){
      
      sidebarMenu(
        menuItem(paste0("Project: ",projectName()), tabName = "projectStatus",
                 icon = icon("chalkboard-teacher"),selected = T),
        lapply(getAnalysisSteps(pipelines[[analysisType()]]),function(x){
          menuItem(x$stepName, tabName = paste0("step_",names(x$stepName)), icon = x$icon)
        }),
        menuItem("Static Report", tabName = "staticProjectReport", icon = icon("file-pdf")),
        menuItem("Dynamic Report", tabName = "dynamicProjectReport", icon = icon("html5"))
        
      )
    }
  })
  
}