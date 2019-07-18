statusbTab<-function(stepTabName){
  tabItem(tabName = stepTabName,
         fluidRow(column(width = 4,
                         gradientBox(title = "Reads before QC", width = 12, icon = 'fa fa-stream',
                                     gradientColor = "blue",  boxToolSize = "xs",
                                     footer = fluidRow(
                                          column(width = 12,
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
                                                               "References: 5")
                                                )
                                          )
                                     
                                         )
                                     )
                             )
         
         
         )
  )
  
}

qcTab<-function(stepTabName){
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