source("Modules/startProject_ui.R")
source("Modules/openProject_ui.R")
source("Modules/changeProjectDir_ui.R")
source("Modules/news_ui.R")

callTabContent<-function(id){
  ns<-NS(id)
  uiOutput(ns("tabContentUI"))
}

tabContentModule<-function(input, output, session, parentSession) {
  ns <- session$ns
  tabButtonLinks<-reactive("")
  
  output$tabContentUI <- renderUI({
    
    if(projectName()!=""){

      pconf<-read.table("projects_eConf.tsv",sep = "\t",header = T,row.names = 1,stringsAsFactors = F)
      tabitems<-list()
      tabitems[["welcome"]]<-tabItem( #fist tab - static
            tabName = "Welcome",
            fluidRow(
              analysisAvailUIm("analysisavails")
            ),
            fluidRow(
              newProjectUIm("startProjectmodule"),
              openProjectUIm("openProjectmodule"),
              changeProjectDirUIm("changeProjectDirmodule")
            ),
            newsUIm("news")
          )
      tabitems[["pstatus"]]<-tabItem(tabName = "projectStatus",
          class = "active",
          fluidRow(
            column(width = 6,
                   gradientBox(title = "Status project",width = 12, icon = 'fa fa-eye',
                               gradientColor = "blue", boxToolSize = "xs", 
                               paste0("Progress of the project according to ",getAnalysisName(pipelines[[analysisType()]]), " pipeline."),
                               footer = fluidRow(column(width = 12,
                                                        progressBar(id = ns("pb_pstatus"), value = as.numeric(pconf["pPercent",]), total = 100, status = "info", 
                                                                    display_pct = TRUE, striped = TRUE, title = projectName())
                               ))# server
                               #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                   )
            ),
            column(width = 6, gradientBox( title = "Last process performed", width = 12,
                                           icon = 'fa fa-tasks', gradientColor = "blue", 
                                           boxToolSize = "xs",
                                           footer = fluidRow(
                                                    if(pconf["lastStep",]=="none"){
                                                      column(width=12,
                                                        tags$h5("Project not started yet"),
                                                        actionBttn(
                                                          inputId = ns("startPipelinebtn"),
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
                                             )
                                           )
                                           # server
                                           #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                )
          ),
          fluidRow(
            lapply(getAnalysisSteps(pipelines[[analysisType()]]),function(x){
                    if(!file.exists(paste0(x$folder,"/",x$stepID,".conf"))){
                      dashblabel<-"not-performed"
                      stepconf<-data.frame()
                      oufolder<-""
                    }else{
                      stepconf<-read.table(paste0(x$folder,"/",x$stepID,".conf"),
                                           sep = "\t", header = T, row.names = 1,
                                           stringsAsFactors = F)
                      dashlabel<-stepconf["statusStep"]
                      oufolder<-paste0(projectName(),"/",x$folder)
                    }
                      gradientBox(title = x$stepName, width = 4, icon = x$iconHTML,
                                  gradientColor = "red",  boxToolSize = "xs",
                                  footer = fluidRow(column(width = 12,
                                      column(width = 12, 
                                             fluidRow("Status: ",getDashboardLabel(dashblabel),
                                                      tags$br(),
                                                      "input files: ",stepconf["niFiles",], tags$br(),
                                                      "output files: ",stepconf["noFiles",], tags$br(),
                                                      "output folder: ",oufolder, tags$br(),
                                                      "Time elapsed: "),stepconf["timeElpased",],
                                             tags$br(),tags$br(),
                                               actionBttn( inputId = ns(paste0(x$stepID,"_tabBtn")),
                                                           label = "Dive into", style = "jelly",
                                                           color = "primary", icon = icon("search")
                                               )

                                      )
                                  ))
                      )     
            })
          )
        )
      ######################## Independent step tabs  ############################################
      for(x in getAnalysisSteps(pipelines[[analysisType()]])){
        #tabitems[[paste0("step_",x$stepID)]]<-x$tabcontent
      }
      tabitems[["step_rbqc"]]<-statusbUIm("statusbmodule","step_rbqc")
      
      ########################################################################################
      tabitems[["preport"]]<-tabItem(tabName = "projectReport",
                fluidRow(column(width = 4,
                     gradientBox(title = "Report summary", width = 12, icon = 'fa fa-file-invoice',
                        gradientColor = "blue",  boxToolSize = "xs",
                        footer = fluidRow(column(width = 12,
                            column(width = 6, 
                                   fluidRow(tags$strong("Static Report")),
                                   fluidRow("Status: done",tags$br(),
                                            "Pages: 10",tags$br(),
                                            "Figures: 3",tags$br(),
                                            "Tables: 3",tags$br(),
                                            "References: 5")
                            ),
                            column(width = 6,
                                   fluidRow(tags$strong("Dynamic Report")),
                                   fluidRow("Status: done",tags$br(),
                                            "Pages: 10",tags$br(),
                                            "Figures: 5",tags$br(),
                                            "Tables: 2",tags$br(),
                                            "References: 5"))
                          )
                        )
                        # server
                        #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                     )
                ),
                column(width=4,
                  gradientBox(title = "Static report options", width = 12, icon = 'fa fa-file-pdf',
                    gradientColor = "blue",  boxToolSize = "xs",
                    footer = fluidRow(
                      column(width = 6,
                          actionBttn( inputId = ns("downStaticRepBtn"),
                                      label = "Save", style = "jelly",
                                      color = "primary", icon = icon("download")
                          )),
                      column(width = 6,
                        actionBttn( inputId = ns("viewStaticRepBtn"),
                                    label = "View", style = "jelly",
                                    color = "primary", icon = icon("eye")
                        )
                      )
                    )
                    # server
                    #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                  )
                ),
                column(width = 4,
                       gradientBox(title = "Dynamic report options", width = 12, icon = 'fa fa-html5',
                           gradientColor = "blue",  boxToolSize = "xs",
                           footer = fluidRow(
                             column(width = 6,
                               actionBttn( inputId = ns("downDynaRepBtn"),
                                           label = "Save", style = "jelly",
                                           color = "primary", icon = icon("download")
                               )),
                             column(width = 6,
                               actionBttn( inputId = ns("viewDynaRepBtn"),
                                           label = "View", style = "jelly",
                                           color = "primary", icon = icon("eye")
                               )
                             )
                           )
                           # server
                           #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                       )
                  )
                ),
                fluidRow(
                  gradientBox(title = "Report view", width = 12, icon = 'fa fa-eye',
                    gradientColor = "blue", boxToolSize = "xs", 
                    footer = fluidRow(column(width = 12

                    )) # server
                    #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                  )
                )
        )
      
      tabitems<-unname(tabitems)
      do.call(tabItems, tabitems)
      
    }else{
      tabItem( #fist tab - static
        tabName = "Welcome",
        fluidRow(
          analysisAvailUIm("analysisavails")
        ),
        fluidRow(
          newProjectUIm("startProjectmodule"),
          openProjectUIm("openProjectmodule"),
          changeProjectDirUIm("changeProjectDirmodule")
        ),
        newsUIm("news")
      )
    }
  })
  

observe({
  if(projectName()=="")return()
  
  lapply(getAnalysisSteps(pipelines[[analysisType()]]), function(x){
    observeEvent(input[[paste0(x$stepID,"_tabBtn")]],{
      
        updateTabItems(session = parentSession, inputId =  "mainMenu",selected =  paste0("step_",x$stepID))
        
      },ignoreNULL = T,ignoreInit = T)

  })
  
  #lapply(getAnalysisSteps(pipelines[[analysisType()]]),function(x){
    callModule(statusbTabModule,"statusbmodule","0-raw","rbqc")
  #})
})





}