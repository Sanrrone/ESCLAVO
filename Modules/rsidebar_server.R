rsidebarModule<-function(input, output, session) {

  output$rsidebarcontent <- renderUI({
    accordion(
      lapply(names(analysisAvail),function(module){
        getAnalysisOptions(pipelines[[module]])
      })
    )
  })
}