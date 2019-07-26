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
                          gradientBox(title = "Software details", width = 12, icon = 'fa fa-barcode',
                                      gradientColor = "purple",  boxToolSize = "xs",
                                      footer = fluidRow(column(width = 12,style = "overflow-x:scroll;",
                                                    tags$strong(soft),sversion,
                                                    tags$br(),
                                                    tags$strong("Work folder: "),folder,
                                                    tags$br(),
                                                    tags$strong("Location: "),spath,
                                                    tags$br()
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
    ffolder<-projectConf()["ffolder",]
    pattern<-projectConf()["fqpattern",]
    if(!file.exists(paste0(ffolder,"/",stepID,".conf"))){
      stepStatus<-"not-performed"
      stepconf<-data.frame(inputFiles=list.files(ffolder,pattern = pattern),
                           qcFile=NA,
                           stringsAsFactors = F)
    }else{
      stepconf<-read.table(paste0(ffolder,"/",stepID,".conf"),
                           sep = "\t", header = T,
                           stringsAsFactors = F)
      stepStatus<-unique(stepconf[,"stepStatus"])
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
      }else if (stepStatus=="done"){
        getDashboardLabel(stepStatus)
      }else{
        column(width=12,
               tags$h4(projectConf()["lastStep",])
        )
      }
    })
    
    output$rbqcTable<-renderDataTable({
      datatable(stepconf,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })
    
    output$readReportSampleUI<-renderUI({
      nsamples<-length(stepconf$qcFile)
      samplesInFolder<-length(stepconf$inputFiles)
      if(nsamples==0 | samplesInFolder == 0){
        gradientBox(title = "Not read status found"," read status", width = 12, icon = 'fa fa-stream',
                    gradientColor = "purple",  boxToolSize = "xs",collapsible = T,
                    footer = fluidRow(column(width=12,tags$h5(paste0("No results are in ",folder))))
        )
      }else{
        addResourcePath("ffolder",ffolder) #add fastq folder to allow iframe load the local html file
        box(title = "Read status report", width = 12, icon = 'fa fa-stream',
            gradientColor = "purple",  boxToolSize = "xs", collapsible = T,
            footer = fluidRow(column(width = 12,
              tags$iframe(seamless = "seamless",
                          src=paste0("ffolder/multiqc_report.html"),
                          height=600, width="100%",frameborder=0)
            ))
        )

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