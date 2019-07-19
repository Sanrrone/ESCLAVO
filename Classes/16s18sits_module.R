statusbUIm<-function(id, stepTabName){
  ns<-NS(id)

  tabItem(tabName = stepTabName,
          fluidRow(column(width = 6,
                          gradientBox(title = "Input files", width = 12, icon = 'fa fa-stream',
                                      gradientColor = "purple",  boxToolSize = "xs",
                                      footer = dataTableOutput(ns("rbqcTable"))
                                      
                          )
                    ),
                    column(width = 6,
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
          uiOutput(ns("readReportSampleUI"))
          
  )
  
}

statusbTabModule<-function(input,output,session, folder, stepID){
  ns<-session$ns
  
  observe({
    if(projectName()=="")return()
    if(!file.exists(paste0(folder,"/",stepID,".conf"))){
      stepStatus<-"not-performed"
      stepconf<-data.frame()
      oufolder<-""
    }else{
      stepconf<-read.table(paste0(folder,"/",stepID,".conf"),
                           sep = "\t", header = T, row.names = 1,
                           stringsAsFactors = F)
      stepStatus<-stepconf["statusStep"]
      oufolder<-paste0(projectName(),"/",folder)
      print(paste0("variables: ",folder," and ",stepID, " and ", oufolder))
      
    }
    
    output$statusbUI<-renderUI({
      print(paste0("variables: ",folder," and ",stepID))
      
      if(stepStatus=="not-performed"){
        column(width=12,
               tags$h5("Reads status step not started yet"),
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
      nsamples<-nrow(stepconf)
      if(nsamples==0)return()
      lapply(1:nsamples, function(i){
        gradientBox(title = paste0(stepconf[i,"sample"]," read status"), width = 12, icon = 'fa fa-stream',
                    gradientColor = "purple",  boxToolSize = "xs",collapsible = T,
                    footer = fluidRow(
                      tags$iframe(src="http://google.cl", height=600, width="100%")
                    )
        )
      })
      
    })

  })
}