tabContentModule<-function(input, output, session) {
  ns <- session$ns
  
  output$tabContentUI <- renderUI({
    if(projectName()!=""){
      tabItems(
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
          ),
        tabItem(tabName = "projectStatus",
          class = "active",
          fluidRow(
            column(width = 6,
              gradientBox(title = "Status project",width = 12, icon = 'fa fa-eye',
                gradientColor = "blue", boxToolSize = "xs", 
                paste0("Progress of the project according to ",getAnalysisName(pipelines[[analysisType()]]), " pipeline."),
                footer = fluidRow(column(width = 12,
                             progressBar(id = "pb1", value = 20, total = 100, status = "info", 
                                         display_pct = TRUE, striped = TRUE, title = projectName())
                      )) # server
                #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
              )
            ),
            column(width = 6, gradientBox( title = "Last process performed", width = 12,
                     icon = 'fa fa-tasks', gradientColor = "blue", 
                     boxToolSize = "xs", "Progress of lastet modified projects",
                     footer = fluidRow(
                       column(width=12,tags$h4("testprocess"))
                     )
                     # server
                     #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                   )
                )
          ),
          fluidRow(
            lapply(getAnalysisSteps(pipelines[[analysisType()]]),function(x){
                      gradientBox(title = "Report summary", width = 4, icon = 'fa fa-file-invoice',
                                  gradientColor = "gray",  boxToolSize = "xs",
                                  footer = fluidRow(column(width = 12,
                                      column(width = 12, 
                                             fluidRow(tags$strong(paste0(x$stepName," summary"))),
                                             fluidRow("Status: ",dashboardLabel("Done", status = "success"),
                                                      tags$br(),
                                                      "input files: 6",tags$br(),
                                                      "output files: 5",tags$br(),
                                                      "output folder: projects/somefolder",tags$br(),
                                                      "Time elapsed: 1 minute"),
                                             tags$br(),tags$br(),
                                               actionBttn( inputId = ns(paste0(x$stepID,"_tabBtn")),
                                                           label = "Dive into", style = "jelly",
                                                           color = "primary", icon = icon("chart-pie")
                                               )

                                      )
                                  )
                                  
                                  )
                      )     
            })
          )
        ),
        tabItem(tabName = "projectReport",
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
      )
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
  

}