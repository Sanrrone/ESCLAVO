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
      stepconf<-read.table(paste0(projectConf()["ffolder",],"/",stepID,".conf"),
                           sep = "\t", header = T,
                           stringsAsFactors = F)
      stepStatus<-unique(stepconf[,"stepStatus"])
      
      output$readReportSamplebUI<-renderUI({
          addResourcePath("ffolder",ffolder) #add fastq folder to allow iframe load the local html file
          box(title = "Read status report", width = 12, icon = 'fa fa-stream',
              gradientColor = "purple",  boxToolSize = "xs", collapsible = T,
              footer = fluidRow(column(width = 12,
                                       tags$iframe(seamless = "seamless",
                                                   src="ffolder/multiqc_report.html",
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
                 inputId = ns("startRSBbtn"),
                 label = "Start",
                 style = "jelly",
                 color = "danger",
                 icon = icon("rocket")
               ))
      }else if (stepStatus=="running"){
        getDashboardLabel(stepStatus)
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
  
  
  observeEvent(input$startRSBbtn,{
    setwd(esclavoHome)
    pfolder<-projectConf()["pfolder",]
    folder<-projectConf()["ffolder",]
    fqpattern<-projectConf()["fqpattern",]
    analysis<-projectConf()["analysis",]#same name for files
    system(paste0("setsid bash pipelines/",analysis,"/",analysis,
                  ".sh --force -p ",pfolder," -f ",folder," -pt ",fqpattern," -m all"),
           wait = F,intern = F, timeout = 0)
  },ignoreNULL = T, ignoreInit = T)
  
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
                 gradientBox(title = "QC step", width = 12, icon = 'fa fa-stream',
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
    pfolder<-paste0(projectConf()["pfolder",],"/1-qc/")
    if(!file.exists(paste0(pfolder,stepID,".conf"))){
      stepStatus<-"not-performed"
      stepconf<-data.frame(inputFiles=list.files(ffolder,pattern = pattern),
                           qcFile=NA,
                           stringsAsFactors = F)
    }else{
      stepconf<-read.table(paste0(pfolder,stepID,".conf"),
                           sep = "\t", header = T,
                           stringsAsFactors = F)
      stepStatus<-unique(stepconf[,"stepStatus"])
      
      output$QCresultsUI<- renderUI({
        if(ncol(stepconf)<=2)return()
        
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
      if(!file.exists(paste0(pfolder,stepID,".conf"))){
        column(width=12,
               tags$h4("Waiting for 'status of reads' step")
        )
      }else{
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
        }else if (stepStatus=="running"){
          getDashboardLabel(stepStatus)
        }else if (stepStatus=="done"){
          getDashboardLabel(stepStatus)
        }else{
          column(width=12,
                 tags$h4(projectConf()["lastStep",])
          )
        }
      }


    })
    
    output$qcReportUI<- renderDataTable({
      if(!file.exists(paste0(pfolder,"qc_summary.tsv")))return()
      qcdf<-read.table(paste0(pfolder,"qc_summary.tsv"),sep = "\t",stringsAsFactors = F)
      datatable(qcdf,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })

  })
  
  observeEvent(input$startQCbtn,{
    setwd(esclavoHome)
    pfolder<-projectConf()["pfolder",]
    folder<-projectConf()["ffolder",]
    fqpattern<-projectConf()["fqpattern",]
    analysis<-projectConf()["analysis",]#same name for files
    system(paste0("setsid bash pipelines/",analysis,"/",analysis,
                  ".sh --force -p ",pfolder," -f ",folder," -pt ",fqpattern," -m 'qc statusa assignTaxonomy report'"),
           wait = F,intern = F, timeout = 0)
  },ignoreNULL = T, ignoreInit = T)
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
    pfolder<-paste0(projectConf()["pfolder",],"/1-qc/")
    
    if(!file.exists(paste0(pfolder,stepID,".conf"))){
      stepStatus<-"not-performed"
      stepconf<-data.frame(inputFiles=list.files(ffolder,pattern = pattern),
                           qcFile=NA,
                           stringsAsFactors = F)
    }else{
      stepconf<-read.table(paste0(pfolder,stepID,".conf"),
                           sep = "\t", header = T,
                           stringsAsFactors = F)
      stepStatus<-unique(stepconf[,"stepStatus"])
      
      output$readReportSampleAfterUI<-renderUI({
        if(ncol(stepconf)<=2 )return()
        
          addResourcePath("ffolder",pfolder) #add fastq folder to allow iframe load the local html file
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
      if(!file.exists(paste0(pfolder,stepID,".conf"))){
        column(width=12,
               tags$h4("Waiting for 'QC' step")
        )
      }else{
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
        }else if (stepStatus=="running"){
          getDashboardLabel(stepStatus)
        }else if (stepStatus=="done"){
          getDashboardLabel(stepStatus)
        }else{
          column(width=12,
                 tags$h4(projectConf()["lastStep",])
          )
        }
      }
    })
    
    output$raqcTable<-renderDataTable({
      datatable(stepconf,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })
    
  })
  
  observeEvent(input$startRSafterbtn,{
    setwd(esclavoHome)
    pfolder<-projectConf()["pfolder",]
    folder<-projectConf()["ffolder",]
    fqpattern<-projectConf()["fqpattern",]
    analysis<-projectConf()["analysis",]#same name for files
    system(paste0("setsid bash pipelines/",analysis,"/",analysis,
                  ".sh --force -p ",pfolder," -f ",folder," -pt ",fqpattern," -m 'statusa assignTaxonomy report'"),
           wait = F,intern = F, timeout = 0)
  },ignoreNULL = T, ignoreInit = T)
}

taxCountTabUIm<-function(id, stepTabName,folder, soft, sversion, spath){
  ns<-NS(id)
  
  tabItem(tabName = stepTabName,
          fluidRow(
            column(width = 6,
                   fluidRow(
                     column(width = 6,
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
                     column(width = 6,
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
                   fluidRow(
                     uiOutput(ns("TCpcaUI"))
                   )
            ),
            column(width = 6,
                   gradientBox(title = "Input files", width = 12, icon = 'fa fa-stream',
                               gradientColor = "purple",  boxToolSize = "xs",
                               footer = dataTableOutput(ns("tcTable"))
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
    pfolder<-paste0(projectConf()["pfolder",],"/2-taxInsight/")
    
    if(!file.exists(paste0(pfolder,stepID,".conf"))){
      stepStatus<-"not-performed"
      stepconf<-data.frame(inputFiles=list.files(ffolder,pattern = pattern),
                           qcFile=NA,
                           stringsAsFactors = F)
    }else{
      stepconf<-read.table(paste0(pfolder,stepID,".conf"),
                           sep = "\t", header = T,
                           stringsAsFactors = F)
      stepStatus<-unique(stepconf[,"stepStatus"])
      if(!file.exists(paste0(pfolder,"tax_table.tsv")))return()
      
      tax<-read.table(paste0(pfolder,"tax_table.tsv"), sep = "\t", header = T,stringsAsFactors = F)
      otu<-read.table(paste0(pfolder,"otu_table.tsv"), sep = "\t", header = T,stringsAsFactors = F)
      tax[is.na(tax)]<-"Unclassified"
      ps<-phyloseq(otu_table(otu, taxa_are_rows=FALSE), tax_table(as.matrix(tax)))
      if(file.exists(paste0(ffolder,'/metadata.tsv'))){
        metadata<-read.table(paste0(ffolder,'/metadata.tsv'), sep='\t', header=T, stringsAsFactors = F)
        rownames(metadata)<-metadata$sample
        sample_data(ps)<-sample_data(metadata)
      }else{
        metadata<-data.frame()
      }
      
      top10 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:10]
      ps.top10 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
      ps.top10 <- prune_taxa(top10, ps.top10)
      
      pst10<-reactiveVal(ps.top10)
      ps<-reactiveVal(ps)
      
      output$TCpcoa <- renderPlotly({
        ordu = ordinate(ps(), 'PCoA', weighted=TRUE)
        plot_ordination(ps(), ordu, color='treatment') + geom_point(size=3)
      })
      
      output$TCresultsUI <- renderUI({
        if(ncol(stepconf)<=2 | !file.exists(paste0(pfolder,"abundance.tsv")))return()
        
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
    
    output$TCpcaUI <- renderUI({
      if(stepStatus!="done")return()
      gradientBox(title = "Principal Component Analysis", width = 12, icon = 'fa fa-stream',
                  gradientColor = "purple",  boxToolSize = "xs",
                  footer = fluidRow(
                    column(width = 12,
                           plotlyOutput(ns("TCpcoa"))
                    )
                  )
      )
    })
    
    
    output$tcUI<-renderUI({
      if(!file.exists(paste0(pfolder,stepID,".conf"))){
        column(width=12,
               tags$h4("Waiting for 'Status after QC' step")
        )
      }else{
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
        }else if (stepStatus=="running"){
          getDashboardLabel(stepStatus)
        }else if (stepStatus=="done"){
          getDashboardLabel(stepStatus)
        }else{
          column(width=12,
                 tags$h4(projectConf()["lastStep",])
          )
        }
      }
    })
    
    output$tcTable<-renderDataTable({
      datatable(stepconf,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })
    
    output$tcabundanceUI<- renderDataTable({
      tcdf<-read.table(paste0(pfolder,"abundance.tsv"),sep = "\t",stringsAsFactors = F,header = T)
      datatable(tcdf,options = list(scrollX=TRUE, scrollCollapse=TRUE))
    })
    
    output$abundanceFamUI<- renderPlot({
      if(nrow(metadata)!=0 & 'treatment' %in% colnames(metadata) & 'sample' %in% colnames(metadata)){
        plot_bar(pst10(), x='Sample', fill='Family') + geom_bar(position='fill', stat='identity', color='black') + 
          theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1)) + facet_wrap(.~treatment,scales = 'free_x')
      }else{
        plot_bar(pst10(), x='Sample', fill='Family') + geom_bar(position='fill', stat='identity', color='black') + 
          theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1))
      }
    })
    
    output$abundanceGenUI<- renderPlot({
      if(nrow(metadata)!=0 & 'treatment' %in% colnames(metadata) & 'sample' %in% colnames(metadata)){
        plot_bar(pst10(), x='Sample', fill='Genus') + geom_bar(position='fill', stat='identity', color='black') + 
          theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1)) + facet_wrap(.~treatment,scales = 'free_x')
      }else{
        plot_bar(pst10(), x='Sample', fill='Genus') + geom_bar(position='fill', stat='identity', color='black') + 
          theme_minimal() + theme(axis.text.x = element_text(angle = 65, hjust = 1))
      }
    })
    
  })
  
  observeEvent(input$startTCbtn,{
    setwd(esclavoHome)
    pfolder<-projectConf()["pfolder",]
    folder<-projectConf()["ffolder",]
    fqpattern<-projectConf()["fqpattern",]
    analysis<-projectConf()["analysis",]#same name for files
    system(paste0("setsid bash pipelines/",analysis,"/",analysis,
                  ".sh --force -p ",pfolder," -f ",folder," -pt ",fqpattern," -m 'assignTaxonomy report'"),
           wait = F,intern = F, timeout = 0)
  },ignoreNULL = T, ignoreInit = T)
}