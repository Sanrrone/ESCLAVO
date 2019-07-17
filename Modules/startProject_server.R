newProjectModule<-function(input, output, session) {
  ns <- session$ns
    
  observeEvent(input$newprojectbtn,{
    showModal(modalDialog(size = "l",
                          easyClose = TRUE,
                          gradientBox(title = "New Project", gradientColor = "red",
                                      closable = F, icon = "fa fa-project-diagram",
                                      width = 12,footer = fluidRow(column(width = 12,
                                        directoryInput(inputId = ns('newProjectFolder'), 
                                                       label = 'Project folder (name)',
                                                       value = projectsFolder),
                                        directoryInput(inputId = ns('newProjectFastq'), 
                                                       label = 'Where are the Project fastq files?'),
                                        uiOutput(ns("fastqOptions")),
                                        awesomeRadio(
                                          inputId = ns("analysisType"),
                                          label = "Select analysis to perform", 
                                          choices = unlist(lapply(names(pipelines),function(x){
                                            niceName<-c(getAnalysisName(pipelines[[x]]))
                                            names(x)<-niceName
                                            x
                                          })),
                                          status = "primary"
                                        ),
                                        tags$br(),
                                        fluidRow(
                                          column(width = 2),
                                          column(width = 3,                  
                                                 actionBttn(
                                                   inputId = ns("submitNewProyect"),
                                                   label = "Go!",
                                                   style = "jelly",
                                                   color = "danger",
                                                   icon = icon("rocket")
                                                 )),
                                          column(width = 1),
                                          column(width = 3,
                                                 modalButton(label = "Dismiss",icon = icon("close"))
                                          )
                                        ),
                                        column(width = 3)
                                      ))

                          ),
                          footer = NULL
            )
    )
  })
  
  observeEvent(eventExpr = input$newProjectFolder, {
                 if(input$newProjectFolder == 0) return()
                 path = choose.dir(default = readDirectoryInput(session, 'newProjectFolder'))
                 # update the widget value
                 updateDirectoryInput(session, 'newProjectFolder', value = path)
  })
  
  observeEvent(eventExpr = input$newProjectFastq, {
    if(input$newProjectFastq == 0){return()}
    fpath = choose.dir(default = readDirectoryInput(session, 'newProjectFastq'))
    if(is.na(fpath)){
      return()
    }else{
      nfastq<-length(list.files(fpath,"*.fastq"))
      if(nfastq == 0){
        sendSweetAlert(session,title = "Warning",
                       text = "No fastq files in the folder selected, please select a folder with fastq files inside.",
                       type = "error")
        return()
      }
    }
    # update the widget value
    updateDirectoryInput(session, 'newProjectFastq', value = fpath)
    
    ppath = readDirectoryInput(session, 'newProjectFolder')
    
    fsplitted<-strsplit(fpath,"/")[[1]]
    psplitted<-strsplit(ppath,"/")[[1]]
    fsplitted<-paste(fsplitted[1:length(fsplitted)-1],collapse = "/")
    psplitted<-paste(psplitted[1:length(psplitted)],collapse = "/")
    fastqoptions<-c("Nothing, just use it" = 1,
                    "Make a copy inside project folder" = 2,
                    "Move entire fastq folder to project folder" = 3)
    if(psplitted == fsplitted){
      fastqoptions<-fastqoptions[1]
    }

    output$fastqOptions <- renderUI({
        fluidRow(
          column(width = 2),
          column(width = 10,
            awesomeRadio(inputId = ns("fastqActions"),label = "What should I do with the fastq?",
                        choices = fastqoptions,
                        status =  "primary",inline = F,
                        selected = 1
            )
          )
        )
    })
    
    
    
  })
  
  observeEvent(input$submitNewProyect,{
    Ppath = readDirectoryInput(session, 'newProjectFolder')
    Fpath = readDirectoryInput(session, 'newProjectFastq')
    if(Fpath == ""){
      sendSweetAlert(session,title = "Warning",text = "No fastq folder select",type = "error")
      return()
    }
    loadingState()
    switch(input$fastqActions,
           "1" = {},
           "2" = {system(paste("cp -r", Fpath, paste0(Ppath,"/0-fastq")),wait = T);Fpath<-paste0(Ppath,"/0-fastq")},
           "3" = {system(paste("mv", Fpath, paste0(Ppath,"/0-fastq")), wait = T);Fpath<-paste0(Ppath,"/0-fastq")}
    )
    setwd(Ppath)
    aoptions<-makeProjectDescription(Fpath,Ppath,
                           input$analysisType,
                           getAnalysisVersion(pipelines[[input$analysisType]]),
                           "open")
    pname<-strsplit(Ppath,"/")[[1]]
    pname<-pname[length(pname)]
    write.table(aoptions,paste0(pname,"_eConf.tsv"),row.names = F,quote = F,sep = "\t")
    #callModule(tabDefModule,"tabdefmodule",pname,input$analysisType)
    projectName(pname)
    analysisType(input$analysisType)
    removeModal()

  })

}
