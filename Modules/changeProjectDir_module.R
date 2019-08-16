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
changeProjectDirModule<-function(input, output, session) {
  ns <- session$ns
  
  observeEvent(input$changeprojectdirbtn,{
    showModal(modalDialog(size = "l",
                          easyClose = TRUE,
                          gradientBox(title = "Choose Projects folder", gradientColor = "purple",
                                      closable = F, icon = "fa fa-open",
                                      width = 12,footer = fluidRow(column(width = 12,
                                                                          directoryInput(inputId = ns('projectsfolder'), 
                                                                                         label = 'Select a folder where are the projects',
                                                                                         value = AllprojectsFolder()),
                                                                          tags$br(),
                                                                          fluidRow(column(width = 2),
                                                                                   column(width = 3,                  
                                                                                          actionBttn(
                                                                                            inputId = ns("newprojectsfolderbtn"),
                                                                                            label = "Go!",
                                                                                            style = "jelly",
                                                                                            color = "danger",
                                                                                            icon = icon("rocket")
                                                                                          )),
                                                                                   column(width = 1),
                                                                                   column(width = 3,
                                                                                          modalButton(label = "Dismiss",icon = icon("close"))
                                                                                   ),
                                                                                   column(width = 3)
                                                                          )
                                      )
                                      
                                      )),
                          footer = NULL
        )
    )
  },ignoreNULL = T, ignoreInit = T)
  
  observeEvent(input$projectsfolder, {
    if(input$projectsfolder == 0) return()
    path = choose.dir(default = readDirectoryInput(session, 'projectsfolder'))
    # update the widget value
    updateDirectoryInput(session, 'projectsfolder', value = path)
  },ignoreNULL = T)
  
  
  observeEvent(input$newprojectsfolderbtn,{
    AllprojectsFolder(readDirectoryInput(session, 'projectsfolder'))
    
    removeModal()
  },ignoreNULL = T, ignoreInit = T)
  
  
  
}


