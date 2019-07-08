changeProjectDirUIm<-function(id){
  ns<-NS(id)
  
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
      style = "jelly",
      color = "warning",
      icon = icon("gears")
    )#,
    #"Begin a new project selecting a folder containing reads."
  )
}