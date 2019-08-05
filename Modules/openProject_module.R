openProjectUIm<-function(id){
  ns<-NS(id)
  
  gradientBox(
    title = "Search existing project", width = 4, icon = "fa fa-folder",
    gradientColor = "red",  boxToolSize = "xs", 
    footer =    "Open an existing project selecting the configuration file of it.", 
    actionBttn(
      inputId = ns("openprojectbtn"), label = "Open", style = "jelly",
      color = "warning", icon = icon("folder-open")
    )
    #"Begin a new project selecting a folder containing reads."
  )
}

openProjectModule<-function(input, output, session, parentSession) {
  ns<-session$ns
  observeEvent(input$openprojectbtn,{
    showModal(modalDialog(size = "l",
                          easyClose = TRUE,
                          gradientBox(title = "Open Project", gradientColor = "purple",
                                      closable = F, icon = "fa fa-open",
                                      width = 12,footer = fluidRow(column(width = 12,
                                                                          directoryInput(inputId = ns('openPfolder'), 
                                                                                         label = 'Project folder where eConf file is',
                                                                                         value = projectsFolder()),
                                                                          tags$br(),
                                                                          fluidRow(column(width = 2),
                                                                            column(width = 3,                  
                                                                                   actionBttn(
                                                                                     inputId = ns("openProjectbtn"),
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
  
  observeEvent(eventExpr = input$openPfolder, {
    if(input$openPfolder == 0) return()
    path = choose.dir(default = readDirectoryInput(session, 'openPfolder'))
    # update the widget value
    updateDirectoryInput(session, 'openPfolder', value = path)
  },ignoreNULL = T)
  
  
  observeEvent(input$openProjectbtn,{
    Ppath<- readDirectoryInput(session, 'openPfolder')
    pname<- strsplit(Ppath,"/")[[1]]
    pname<- pname[length(pname)]
    if(!file.exists(paste0(Ppath,"/",paste0(pname,"_eConf.tsv")))){
      sendSweetAlert(session,title = "Warning",text = "No ESCLAVO config file found",type = "error")
      return()
    }
    setwd(Ppath)
    aoptions<-read.csv(paste0(pname,"_eConf.tsv"),header = T,sep = "\t",stringsAsFactors = F)
    projectConf(aoptions)
    analysisType(projectConf()["analysis",])
    projectFolder(Ppath)
    projectName(pname)
    
    removeModal()
  },ignoreNULL = T, ignoreInit = T)
  
}