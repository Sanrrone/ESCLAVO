analysisAvailUIm <- function(id,analysisAvail) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  uiOutput(ns("analysisImages"))
  
}

analysisAvailsModule<-function(input, output, session) {

  output$analysisImages<-renderUI({
    widgetUserBox(
      title = "Available analysis",
      subtitle = "Unlocked",
      width = 12,
      type = 2,
      src = dataURI(file="images/unlocked.ico", mime="image/png"),
      color = "aqua-active",
      userList(
        lapply(names(analysisAvail),function(module){
          getAnalysisImages(pipelines[[module]])
        })
      )
    )

  })
}