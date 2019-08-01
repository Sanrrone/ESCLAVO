source("Modules/changeProjectDir_ui.R")
source("Modules/news_ui.R")

callTabContent<-function(id){
  ns<-NS(id)
  uiOutput(ns("tabContentUI"))
}

tabContentModule<-function(input, output, session, parentSession) {
  
  ns <- session$ns
  tabButtonLinks<-reactive("")

  observe({
    if(nrow(projectConf())==0){
      output$tabContentUI <- renderUI({
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
     })
    }else{
      atype<-projectConf()["analysis",]
      output$tabContentUI <- renderUI({
        pipelines[[atype]]<-setNewStepFolder(pipelines[[atype]],1,projectConf()["ffolder",])
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
                                                            paste0("Progress of the project according to ",getAnalysisName(pipelines[[atype]]), " pipeline."),
                                                            footer = fluidRow(column(width = 12,
                                                                                     progressBar(id = ns("pb_pstatus"), value = as.numeric(projectConf()["pPercent",]), 
                                                                                                 total = 100, status = "info", 
                                                                                                 display_pct = TRUE, striped = TRUE, title = projectName())
                                                            ))# server
                                                            #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                                                )
                                         ),
                                         column(width = 6, gradientBox( title = "Last process performed", width = 12,
                                                                        icon = 'fa fa-tasks', gradientColor = "blue", 
                                                                        boxToolSize = "xs",
                                                                        footer = fluidRow(
                                                                          if(projectConf()["lastStep",]=="none"){
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
                                                                                   tags$strong("Last step: "),projectConf()["lastStep",]
                                                                            )
                                                                          }
                                                                        )
                                         )
                                         # server
                                         #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                                         )
                                       ),
                                       fluidRow(
                                         column(width=12,
                                                gradientBox(title = "CPU usage",width = 12, icon = 'fa fa-microchip',
                                                            gradientColor = "red", boxToolSize = "xs", 
                                                            footer = fluidRow(column(width = 12,
                                                                                     plotOutput(ns("cpuUsageBox"))
                                                            ))# server
                                                            #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                                                )
                                         )
                                       ),
                                       fluidRow(
                                         lapply(getAnalysisSteps(pipelines[[atype]]),function(x){
                                           confFileName<-paste0(x$folder,"/",x$stepID,".conf")
                                           if(!file.exists(confFileName)){
                                             dashblabel<-"not-performed"
                                             stepconf<-data.frame(inputFiles=list.files(x$folder,projectConf()["fqpattern",]),
                                                                  timeElpased=rep("0:0:0",length(list.files(x$folder,projectConf()["fqpattern",])))
                                                                  ,stringsAsFactors = F)
                                             
                                           }else{
                                             stepconf<-read.table(confFileName, sep = "\t", header = T, stringsAsFactors = F)
                                             dashblabel<-unique(stepconf[,"stepStatus"])
                                             
                                           }
                                           gradientBox(title = x$stepName, width = 4, icon = x$iconHTML,
                                                       gradientColor = "red",  boxToolSize = "xs",
                                                       footer = fluidRow(column(width = 12,
                                                                                column(width = 12, 
                                                                                       fluidRow("Status: ",getDashboardLabel(dashblabel),
                                                                                                tags$br(),
                                                                                                "input files: ",nrow(stepconf), tags$br(),
                                                                                                "output files: ",length(list.files(x$folder))-nrow(stepconf), tags$br(),
                                                                                                "output folder: ",x$folder, tags$br(),
                                                                                                "Time elapsed: ",parseTimes(stepconf[,"timeElpased"],"minutes")),
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
        for(x in getAnalysisSteps(pipelines[[atype]])[1]){
          tabitems[[paste0("step_",x$stepID)]]<-x$tabcontentUI$ui(x$tabcontentUI$id,
                                                                  paste0("step_",x$stepID), projectConf()["ffolder",],
                                                                  x$software, x$version, x$spath)
        }
        
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
        
      })
      
      
      lapply(getAnalysisSteps(pipelines[[atype]]), function(x){
        observeEvent(input[[paste0(x$stepID,"_tabBtn")]],{
          
          updateTabItems(session = parentSession, inputId =  "mainMenu", selected =  paste0("step_",x$stepID))
          
        },ignoreNULL = T,ignoreInit = T)
        
      })
      
      
      
      cputop<-reactivePoll(intervalMillis = 2000, session = session,
                           checkFunc = function(){ digest("/proc/stat",algo="md5",file=TRUE) },
                           valueFunc = function(){ getCPUusage() })
      
      output$cpuUsageBox <- renderCachedPlot(expr = {
        #options(warn=-1)
        ggplot(cputop(),aes(x=cpu,y=usage,fill=cpu)) + geom_bar(stat = "identity") +
          ylim(0,110) + theme_minimal() +
          geom_text(data = cputop(), nudge_y = 10, angle = 0,
                    aes(x = cpu, y = usage, label = usage))
        #options(warn=0)
      },cputop())
      

    }
    
  })
  
  #check when start pipeline is pressed
  observeEvent(input$startPipelinebtn,{
    setwd(esclavoHome)
    pfolder<-projectConf()["pfolder",]
    folder<-projectConf()["ffolder",]
    fqpattern<-projectConf()["fqpattern",]
    analysis<-projectConf()["analysis",]#same name for files
    removeUI(selector='#startPipelinebtn', immediate=TRUE)
    system2(command = "bash",
            args = paste0("pipelines/",analysis,".sh --force -p ",pfolder," -f ",folder," -pt ",fqpattern),
            wait = F)
    removeUI(selector='#startPipelinebtn', immediate=TRUE)
    
    
  },ignoreNULL = T, ignoreInit = T)

}