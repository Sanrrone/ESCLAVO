openProjectUIm<-function(id){
  ns<-NS(id)
  
  gradientBox(
    title = "Search existing project",
    width = 4,
    icon = "fa fa-folder",
    gradientColor = "red", 
    boxToolSize = "xs", 
    footer =    "Open an existing project selecting the configuration file of it.", 
    actionBttn(
      inputId = ns("openprojectbtn"),
      label = "Open",
      style = "jelly",
      color = "warning",
      icon = icon("folder-open")
    )#,
    #"Begin a new project selecting a folder containing reads."
    
  )
}