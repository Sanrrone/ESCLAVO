newProjectUIm<-function(id){
  ns<-NS(id)
  
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
      style = "jelly",
      color = "warning",
      icon = icon("rocket")
    )#,
    #"Begin a new project selecting a folder containing reads."
  )
}