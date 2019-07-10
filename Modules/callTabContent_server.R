tabContentModule<-function(input, output, session) {
  ns <- session$ns
  
  output$tabContentUI <- renderUI({
    if(projectName()!=""){
      
      tabItem(tabName = ns("projectStatus"),class = "active",
        fluidRow(
          column(width = 6,
            gradientBox(
              title = "Status progress",
              width = 12,
              icon = 'fa fa-tasks',
              gradientColor = "blue", 
              boxToolSize = "xs", 
              "Progress of lastet modified projects",
              footer = fluidRow(column(width = 12,
                           progressBar(id = "pb1", value = 20, total = 100, status = "info", 
                                       display_pct = TRUE, striped = TRUE, title = projectName()),
                    )) # server
              #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
            )
          ),
          column(width = 6, gradientBox(
                   title = "Last process performed",
                   width = 12,
                   icon = 'fa fa-tasks',
                   gradientColor = "blue", 
                   boxToolSize = "xs", 
                   "Progress of lastet modified projects",
                   footer = fluidRow(
                     tags$h4("testprocess")
                   )
                   # server
                   #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                 )
              )
        ),
        fluidRow(
          
        )
      )
      #lapply(getAnalysisSteps(pipelines[[analysisType()]]),function(x){
      #  tabItem(tabName = paste0("step_",names(x$stepName))
      #  )
      #})
      #tabItem(tabName = "staticProjectReport"
      #)
      #tabItem(tabName = "dynamicProjectReport"
      #)
    }
  })
  

}