

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