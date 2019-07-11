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
              gradientBox(title = "Status progress",width = 12, icon = 'fa fa-tasks',
                gradientColor = "blue", boxToolSize = "xs", 
                "Progress of lastet modified projects",
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
          )
        ),
        #lapply(getAnalysisSteps(pipelines[[analysisType()]]),function(x){
        #  tabItem(tabName = paste0("step_",names(x$stepName))
        #  )
        #}),
        tabItem(tabName = "staticProjectReport",
                fluidRow(
                  gradientBox(title = "Static report options", width = 12, icon = 'fa file-pdf',
                    gradientColor = "red",  boxToolSize = "xs", 
                    footer = fluidRow(column(width = 12,
                               actionBttn( inputId = ns("downStaticRepBtn"),
                                 label = "Download report", style = "jelly",
                                 color = "success", icon = icon("download")
                               ),
                               actionBttn( inputId = ns("genStaticRepBtn"),
                                 label = "Re-generate report", style = "jelly",
                                 color = "success", icon = icon("download")
                               )
                    )) # server
                    #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                  )
                ),
                fluidRow(
                  gradientBox( title = "Static report", width = 12, icon = 'fa fa-tasks',
                    gradientColor = "blue", boxToolSize = "xs", 
                    footer = fluidRow(column(width = 12

                    )) # server
                    #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                  )
                )
        ),
        tabItem(tabName = "dynamicProjectReport",
                fluidRow(
                  gradientBox(title = "Dynamic report options", width = 12, icon = 'fa html5',
                              gradientColor = "blue",  boxToolSize = "xs", 
                              footer = fluidRow(column(width = 12,
                                       actionBttn( inputId = ns("downDynamicRepBtn"),
                                                   label = "Download report", style = "jelly",
                                                   color = "success", icon = icon("download")
                                       ),
                                       actionBttn( inputId = ns("genDynamicRepBtn"),
                                                   label = "Re-generate report", style = "jelly",
                                                   color = "success", icon = icon("download")
                                       )
                              )) # server
                              #updateProgressBar(session = session, id = "pb8", value = input$slider, total = 5000)
                  )
                ),
                fluidRow(
                  gradientBox( title = "Dynamic report", width = 12, icon = 'fa fa-tasks',
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