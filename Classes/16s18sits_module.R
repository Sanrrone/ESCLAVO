statusbUIm<-function(id, stepTabName,folder, soft, sversion, spath){
  ns<-NS(id)

  tabItem(tabName = stepTabName,
          fluidRow(column(width = 5,
                          gradientBox(title = "Input files", width = 12, icon = 'fa fa-stream',
                                      gradientColor = "purple",  boxToolSize = "xs",
                                      footer = dataTableOutput(ns("rbqcTable"))
                          )
                    ),
                   column(width = 4,
                          gradientBox(title = "Software details", width = 12, icon = 'fa scroll',
                                      gradientColor = "purple",  boxToolSize = "xs",
                                      footer = fluidRow(column(width = 12,
                                                    tags$strong(soft),sversion,
                                                    tags$br(),
                                                    tags$strong("Work folder: "),folder,
                                                    tags$br(),
                                                    tags$strong("Location: "),spath
                                                  )
                                                )
                                      )
                   ),
                   column(width = 3,
                           gradientBox(title = "Reads before QC", width = 12, icon = 'fa fa-stream',
                                       gradientColor = "purple",  boxToolSize = "xs",
                                       footer = fluidRow(
                                         column(width = 12,
                                            uiOutput(ns("statusbUI"))
                                         )
                                       )
                           )
                    )
          ),
          fluidRow(
            uiOutput(ns("readReportSampleUI"))
          )
          
  )
  
}
statusbTabModule<-function(input,output,session, stepID){
  ns <- session$ns
  
  observe({
    if(projectName()=="")return()
    folder<-projectConf()["ffolder",]
    if(!file.exists(paste0(folder,"/",stepID,".conf"))){
      stepStatus<-"not-performed"
      stepconf<-data.frame(inputFiles=list.files(folder,pattern = ".fastq"),
                           multiqcFile=NA,
                           stringsAsFactors = F)
      oufolder<-""
    }else{
      stepconf<-read.table(paste0(folder,"/",stepID,".conf"),
                           sep = "\t", header = T, row.names = 1,
                           stringsAsFactors = F)
      stepStatus<-stepconf["statusStep"]
      oufolder<-paste0(projectName(),"/",folder)

    }
    
    output$statusbUI<-renderUI({
      
      if(stepStatus=="not-performed"){
        column(width=12,
               tags$h5("Reads status not started yet"),
               actionBttn(
                 inputId = ns("startRSbtn"),
                 label = "Start",
                 style = "jelly",
                 color = "danger",
                 icon = icon("rocket")
               ))
      }else{
        column(width=12,
               tags$h4(pconf["lastStep",])
        )
      }
    })
    
    output$rbqcTable<-renderDataTable({
      return(stepconf)
    })
    
    output$readReportSampleUI<-renderUI({
      nsamples<-length(stepconf$multiqcFile)
      samplesInFolder<-length(stepconf$sample)
      if(nsamples==0 | samplesInFolder == 0){
        gradientBox(title = "Not read status found"," read status", width = 12, icon = 'fa fa-stream',
                    gradientColor = "purple",  boxToolSize = "xs",collapsible = T,
                    footer = fluidRow(column(width=12,tags$h5(paste0("No results are in ",folder))))
        )
      }else{
        
        box(title = paste0(stepconf[1,"sample"]," read status"), width = 12, icon = 'fa fa-stream',
                      gradientColor = "purple",  boxToolSize = "xs",collapsible = T,
                      footer = fluidRow(
                        tags$iframe(src="http://google.cl", height=600, width="100%")
                      )
        )
        if(nsamples>=2){
          lapply(2:nsamples, function(i){
            box(title = paste0(stepconf[i,"sample"]," read status"), width = 12, icon = 'fa fa-stream',
                        status = "info", collapsible = T,collapsed=T,
                        footer = fluidRow(
                          tags$iframe(src="http://google.cl", height=600, width="100%")
                        )
            )
          })
        }

      }

      
    })

  })
}

qcUIm<-function(id,stepTabName, folder, soft, sversion, spath){
  tabItem(tabName = stepTabName,
          fluidRow(column(width = 4,
                          gradientBox(title = "Quality Control", width = 12, icon = 'fa fa-stream',
                                      gradientColor = "blue",  boxToolSize = "xs",
                                      footer = fluidRow(
                                      )
                          )
          )
          
          
          )
  )
  
}
qcModule<-function(input,output,session,folder,stepID)
statusaTab<-function(stepTabName){
  tabItem(tabName = stepTabName,
          fluidRow(column(width = 4,
                          gradientBox(title = "Reads after QC", width = 12, icon = 'fa fa-stream',
                                      gradientColor = "blue",  boxToolSize = "xs",
                                      footer = fluidRow(
                                      )
                          )
          )
          
          
          )
  )
  
}

taxCountTab<-function(stepTabName){
  tabItem(tabName = stepTabName,
          fluidRow(column(width = 4,
                          gradientBox(title = "Taxonomic insight", width = 12, icon = 'fa fa-stream',
                                      gradientColor = "blue",  boxToolSize = "xs",
                                      footer = fluidRow(
                                      )
                          )
          )
          
          
          )
  )
  
}