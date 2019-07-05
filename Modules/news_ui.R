newsUIm<-function(id){
  ns <- NS(id)
  fluidRow(
    gradientBox(
      title = "ESCLAVO news",
      width = 6,
      icon = "fa fa-newspaper",
      gradientColor = "green", 
      boxToolSize = "xs",
      "Check the last analysis updates",
      footer = fluidRow(column(width = 12, 
      timelineBlock(
        timelineEnd(color = "danger"),
        timelineLabel("August 2019", color = "teal"),
        timelineItem(
          title = "First release",
          icon = "rocket",
          color = "olive",
          time = "now",
          "Release with metagenomic and 16S/18S pipelines"
        ),
        timelineItem(
          title = "Centrifuge included",
          footer = "more analysis soon",
          border = FALSE
        ),
        timelineStart(color = "gray")
      )))
    ),
    gradientBox(
      title = "Project progress",
      width = 6,
      icon = 'fa fa-tasks',
      gradientColor = "green", 
      boxToolSize = "xs", 
      "Progress of lastet modified projects",
      footer = fluidRow(column(width = 12,
        progressBar(id = "pb1", value = 20, total = 100, status = "info", display_pct = TRUE, striped = TRUE, title = "Project AG1501"),
        progressBar(id = "pb2", value = 40, total = 100, status = "info", display_pct = TRUE, striped = TRUE, title = "Project AG1820"),
        progressBar(id = "pb2", value = 80, total = 100, status = "info", display_pct = TRUE, striped = TRUE, title = "Project AG0111")
      )
      )
      # ui

      # server
      #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
      
      
    )
  )
  
  
}