newsUIm<-function(id){
  ns <- NS(id)
  
  fluidRow(
    gradientBox(
      title = "ESCLAVO news", width = 6, icon = "fa fa-newspaper", gradientColor = "green",  boxToolSize = "xs",
      "Check the last analysis updates",
      footer = fluidRow(column(width = 12, 
      timelineBlock(
        timelineEnd(color = "danger"), timelineLabel("August 2019", color = "teal"),
        timelineItem(
          title = "First release", icon = "rocket", color = "olive", time = "this month",
          "Release 16S/18S pipeline"
        ),
        timelineItem(
          title = "DADA2 included as Taxonomic classification", footer = "more analysis soon", border = FALSE
        ),
        timelineStart(color = "gray")
      )))
    ),
    gradientBox(
      title = "Project progress", width = 6, icon = 'fa fa-tasks', gradientColor = "green", boxToolSize = "xs", 
      "Progress of lastet modified projects",
      footer = fluidRow(
          uiOutput(ns("projectProgressUI"))
      )
    )
  )
  
}
newsServerm<-function(input, output, session, parentSession) {
  ns <- session$ns
  
  observe({
    if(AllprojectsFolder()==""){
      output$projectProgressUI<- renderUI({
          column(width = 12,
                 progressBar(id = "pb1", value = 20, total = 100, status = "info", display_pct = TRUE, striped = TRUE, title = "Example project 1"),
                 progressBar(id = "pb2", value = 40, total = 100, status = "info", display_pct = TRUE, striped = TRUE, title = "Example project 2"),
                 progressBar(id = "pb2", value = 80, total = 100, status = "info", display_pct = TRUE, striped = TRUE, title = "Example prpject 3")   
          )
      })
    
    }else{
      output$projectProgressUI<- renderUI({
        
        esclavoConfigFiles<-list()
        for(dir in list.files(AllprojectsFolder(),all.files = F,recursive = F)){
          if(file.exists(paste0(dir,"/",dir,"_eConf.tsv"))){
            esclavoConfigFiles[[dir]]<-paste0(dir,"/",dir,"_eConf.tsv")
          }
        }
        
        column(width = 12,
               lapply(esclavoConfigFiles, function(x){
                 
                   aoptions<-read.csv(x, header = T,sep = "\t",stringsAsFactors = F,row.names = 1)
                   pname<-aoptions["pname",]
                   progressBar(id = pname, value = as.numeric(aoptions["pPercent",]), total = 100, status = "info", 
                               display_pct = TRUE, striped = TRUE, title = paste0("Project: ",pname))
                 
               })
        )

      })

    }
    
    


  })
  

}