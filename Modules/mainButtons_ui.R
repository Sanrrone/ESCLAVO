mainButtonUIm<-function(id){
  ns<-NS(id)
  
  fluidRow(
    
    gradientBox(
      title = "Create new project",
      width = 4,
      icon = "fa fa-magic",
      gradientColor = "red", 
      boxToolSize = "xs", 
      footer =    "Begin a new project selecting a folder containing reads.", 
      actionBttn(
        inputId = ns("newprojectbtn"),
        label = "Start",
        style = "material-flat",
        color = "warning",
        icon = icon("rocket")
      )#,
      #"Begin a new project selecting a folder containing reads."

    ),
    gradientBox(
      title = "Search existing project",
      width = 4,
      icon = "fa fa-folder",
      gradientColor = "red", 
      boxToolSize = "xs", 
      footer =    "Open an existing project selecting the Configuration file of it.", 
      actionBttn(
        inputId = ns("openprojectbtn"),
        label = "Open",
        style = "material-flat",
        color = "warning",
        icon = icon("folder-open")
      )#,
      #"Begin a new project selecting a folder containing reads."
      
    ),
    gradientBox(
      title = "Change project directory",
      width = 4,
      icon = "fa fa-people-carry",
      gradientColor = "red", 
      boxToolSize = "xs", 
      footer =    "Change project directory to a folder who contains others project folders", 
      actionBttn(
        inputId = ns("changeprojectdirbtn"),
        label = "Choose",
        style = "material-flat",
        color = "warning",
        icon = icon("gears")
      )#,
      #"Begin a new project selecting a folder containing reads."
      
    )
  )
}
