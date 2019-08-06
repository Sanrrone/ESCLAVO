statusbUIm<-function(id, stepTabName,folder, soft, sversion, spath){
  ns<-NS(id)

  tabItem(tabName = stepTabName,
          fluidRow(column(width = 6,
                          gradientBox(title = "Input files", width = 12, icon = 'fa fa-stream',
                                      gradientColor = "purple",  boxToolSize = "xs",
                                      footer = dataTableOutput(ns("rbqcTable"))
                          )
                    ),
                   column(width = 3,
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
            uiOutput(ns("readReportSamplebUI"))
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
      
      output$readReportSamplebUI<-renderUI({
          addResourcePath("ffolder",ffolder) #add fastq folder to allow iframe load the local html file
          box(title = "Read status report", width = 12, icon = 'fa fa-stream',
              gradientColor = "purple",  boxToolSize = "xs", collapsible = T,
              footer = fluidRow(column(width = 12,
                                       tags$iframe(seamless = "seamless",
                                                   src=paste0("ffolder/multiqc_report.html"),
                                                   height=600, width="100%",frameborder=0)
              ))
          )
      })
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

  })
}

qcUIm<-function(id, stepTabName, folder, soft, sversion, spath){
  ns<-NS(id)
  
  tabItem(tabName = stepTabName,
          fluidRow(column(width = 6,
                          gradientBox(title = "Input files", width = 12, icon = 'fa fa-stream',
                                      gradientColor = "purple",  boxToolSize = "xs",
                                      footer = dataTableOutput(ns("qcTable"))
                          )
          ),
          column(width = 3,
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
                                      uiOutput(ns("qcUI"))
                               )
                             )
                 )
          )
          ),
          uiOutput(ns("QCresultsUI"))

  )
  
}
qcModule<-function(input,output,session,stepID){
  ns <- session$ns
  
  observe({
    if(projectName()=="")return()
    
    ffolder<-projectConf()["ffolder",]
    pattern<-projectConf()["fqpattern",]
    if(!file.exists(paste0("1-qc/",stepID,".conf"))){
      stepStatus<-"not-performed"
      stepconf<-data.frame(inputFiles=list.files(ffolder,pattern = pattern),
                           qcFile=NA,
                           stringsAsFactors = F)
    }else{
      stepconf<-read.table(paste0("1-qc/",stepID,".conf"),
                           sep = "\t", header = T,
                           stringsAsFactors = F)
      stepStatus<-unique(stepconf[,"stepStatus"])
      
      output$QCresultsUI<- renderUI({
        fluidRow(gradientBox(title = "QC Summary", width = 12, icon = 'fa fa-stream',
                             gradientColor = "purple",  boxToolSize = "xs",
                             footer = dataTableOutput(ns("qcReportUI"))
              )
        )
      })
    }
    
    output$qcTable<-renderDataTable({
      datatable(stepconf,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })
    
    output$qcUI<-renderUI({
      
      if(stepStatus=="not-performed"){
        column(width=12,
               tags$h5("QC not started yet"),
               actionBttn(
                 inputId = ns("startQCbtn"),
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
    

    output$qcReportUI<- renderDataTable({
      qcdf<-read.table(paste0("1-qc/qc_summary.tsv"),sep = "\t",stringsAsFactors = F)
      datatable(qcdf,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })

    
  })
}

  
statusaUIm<-function(id, stepTabName,folder, soft, sversion, spath){
  ns<-NS(id)
  
  tabItem(tabName = stepTabName,
          fluidRow(column(width = 6,
                          gradientBox(title = "Input files", width = 12, icon = 'fa fa-stream',
                                      gradientColor = "purple",  boxToolSize = "xs",
                                      footer = dataTableOutput(ns("raqcTable"))
                          )
          ),
          column(width = 3,
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
                 gradientBox(title = "Reads after QC", width = 12, icon = 'fa fa-stream',
                             gradientColor = "purple",  boxToolSize = "xs",
                             footer = fluidRow(
                               column(width = 12,
                                      uiOutput(ns("statusaUI"))
                               )
                             )
                 )
          )
          ),
          fluidRow(
            uiOutput(ns("readReportSampleAfterUI"))
          )
          
  )
  
}
statusaTabModule<-function(input,output,session, stepID){
  ns <- session$ns
  
  observe({
    if(projectName()=="")return()
    ffolder<-projectConf()["ffolder",]
    pattern<-projectConf()["fqpattern",]
    if(!file.exists(paste0("1-qc/",stepID,".conf"))){
      stepStatus<-"not-performed"
      stepconf<-data.frame(inputFiles=list.files(ffolder,pattern = pattern),
                           qcFile=NA,
                           stringsAsFactors = F)
    }else{
      stepconf<-read.table(paste0("1-qc/",stepID,".conf"),
                           sep = "\t", header = T,
                           stringsAsFactors = F)
      stepStatus<-unique(stepconf[,"stepStatus"])
      output$readReportSampleAfterUI<-renderUI({
          addResourcePath("ffolder","1-qc") #add fastq folder to allow iframe load the local html file
          box(title = "Read status report", width = 12, icon = 'fa fa-stream',
              gradientColor = "purple",  boxToolSize = "xs", collapsible = T,
              footer = fluidRow(column(width = 12,
                                       tags$iframe(seamless = "seamless",
                                                   src=paste0("ffolder/multiqc_report.html"),
                                                   height=600, width="100%",frameborder=0)
              ))
          )
      })
      
    }
    
    output$statusaUI<-renderUI({
      
      if(stepStatus=="not-performed"){
        column(width=12,
               tags$h5("Reads status not started yet"),
               actionBttn(
                 inputId = ns("startRSafterbtn"),
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
    
    output$raqcTable<-renderDataTable({
      datatable(stepconf,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })
    
  })
}

taxCountTabUIm<-function(id, stepTabName,folder, soft, sversion, spath){
  ns<-NS(id)
  
  tabItem(tabName = stepTabName,
          fluidRow(column(width = 6,
                          gradientBox(title = "Input files", width = 12, icon = 'fa fa-stream',
                                      gradientColor = "purple",  boxToolSize = "xs",
                                      footer = dataTableOutput(ns("tcTable"))
                          )
          ),
          column(width = 3,
                 gradientBox(title = "Software details", width = 12, icon = 'fa fa-barcode',
                             gradientColor = "purple",  boxToolSize = "xs",
                             footer = fluidRow(column(width = 12, style = "overflow-x:scroll;",
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
                 gradientBox(title = "Taxonomic count step", width = 12, icon = 'fa fa-stream',
                             gradientColor = "purple",  boxToolSize = "xs",
                             footer = fluidRow(
                               column(width = 12,
                                      uiOutput(ns("tcUI"))
                               )
                             )
                 )
          )
          ),
          uiOutput(ns("TCresultsUI"))
          
  )
  
}
taxCountTabModule<-function(input,output,session, stepID){
  ns <- session$ns
  
  observe({
    if(projectName()=="")return()
    ffolder<-projectConf()["ffolder",]
    pattern<-projectConf()["fqpattern",]
    if(!file.exists(paste0("2-taxInsight/",stepID,".conf"))){
      stepStatus<-"not-performed"
      stepconf<-data.frame(inputFiles=list.files(ffolder,pattern = pattern),
                           qcFile=NA,
                           stringsAsFactors = F)
    }else{
      stepconf<-read.table(paste0("2-taxInsight/",stepID,".conf"),
                           sep = "\t", header = T,
                           stringsAsFactors = F)
      stepStatus<-unique(stepconf[,"stepStatus"])
      tax<-read.table("2-taxInsight/tax_table.tsv",sep = "\t",header = T,stringsAsFactors = F)
      otu<-read.table("2-taxInsight/otu_table.tsv",sep = "\t",header = T,stringsAsFactors = F)
      
      ps<-phyloseq(otu_table(otu, taxa_are_rows=FALSE), tax_table(as.matrix(tax)))
      top10 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:10]
      ps.top10 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
      ps.top10 <- prune_taxa(top10, ps.top10)
      ps<-reactiveVal(ps.top10)
      
      output$TCresultsUI <- renderUI({
        fluidRow(
          fluidRow(gradientBox(title = "Taxonomic abundance", width = 12, icon = 'fa fa-calculator',
                               gradientColor = "purple",  boxToolSize = "xs",
                               footer = dataTableOutput(ns("tcabundanceUI"))
          )
          ),
          fluidRow(
            column(width = 6,
                   gradientBox(title = "Top 10 abundance at family level", width = 12, icon = 'fa fa-calculator',
                               gradientColor = "purple",  boxToolSize = "xs",
                               footer = plotOutput(ns("abundanceFamUI"))
                   )
            ),
            column(width = 6,
                   gradientBox(title = "Top 10 abundance at genus level", width = 12, icon = 'fa fa-calculator',
                               gradientColor = "purple",  boxToolSize = "xs",
                               footer = plotOutput(ns("abundanceGenUI"))
                   )
            )
          )
        )

      })

    }
    
    output$tcUI<-renderUI({
      
      if(stepStatus=="not-performed"){
        column(width=12,
               tags$h5("Reads status not started yet"),
               actionBttn(
                 inputId = ns("startTCbtn"),
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
    
    output$tcTable<-renderDataTable({
      datatable(stepconf,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })
    
    output$tcabundanceUI<- renderDataTable({
      tcdf<-read.table(paste0("2-taxInsight/abundance.tsv"),sep = "\t",stringsAsFactors = F,header = T)
      datatable(tcdf,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })
    
    output$abundanceFamUI<- renderPlot({
      plot_bar(ps(), x='Sample', fill='Family') + geom_bar(position='fill', stat='identity', color='black') + 
        theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1))
    })
    
    output$abundanceGenUI<- renderPlot({
      plot_bar(ps(), x='Sample', fill='Genus') + geom_bar(position='fill', stat='identity', color='black') + 
        theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1))
    })
    
  })
}