# Module UI function
headerUIm <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  dashboardHeaderPlus(
    title = tagList(
      span(class = "logo-lg", "ESCLAVO"),icon("heart")),
    enable_rightsidebar = T,
    rightSidebarIcon = "gears",
    left_menu = tagList(
      #uiOutput(ns("header_proyectName"))
      #uiOutput(ns("header_proyectType"))
      #uiChangeThemeDropdown()
      uiOutput("fakefake")
    )
  )
}

headerModule<-function(input, output, session) {
  ns <- session$ns
  
  observe({
    if(projectName()=="")return()
    
    output$header_proyectName <- renderUI({
      dropdownBlock(
        id = ns("headerPname"),
        title = projectName(),
        icon = "sliders",
        selectInput("executeProject",label="Run type",
                    choices = c("Local"=1,"Slurm (not available)"=2),selected = 1)
      )
    })
    
    output$header_proyectType <- renderUI({
      dropdownBlock(
        id = ns("headerPtype"),
        title = analysisType(),
        icon = "sliders"
      )
    })
  })

}