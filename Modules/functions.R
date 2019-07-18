
#analysisAvail=c("16s18sits"=1,"rnaseq"=2,"metagenomic"=3,"gbs"=4)
analysisAvail=c("16s18sits"=1,"metagenomic"=3)

pipelines<-list()
for(module in names(analysisAvail)){
  source(paste0("Classes/",module,"_object.R"))
  pipelines[[module]]<-getpipeline()
}
projectName<-reactiveVal("")
analysisType<-reactiveVal("")
projectsFolder="/home/sandro/Programas/ESCLAVO/projects"

makeProjectDescription<-function(projectFolder,fastqFolder,analysisType, aVersion, pStatus,pPercent,lStep){
  rnames<-c("pfolder","ffolder","analysis","aversion","created","status","pPercent","lastStep")
  
  df<-data.frame(value=c(projectFolder,
                                       fastqFolder,analysisType,aVersion,
                                       Sys.Date(),pStatus,pPercent,lStep),
                 stringsAsFactors = F)
  rownames(df)<-rnames
  
  
  return(df)
}

getDashboardLabel<-function(stepStatus){
  switch(stepStatus,
         "not-performed" = dashboardLabel("Not performed", status = "warning"),
         "running" = dashboardLabel("Running", status = "info"),
         "error" = dashboardLabel("Error", status = "danger"),
         "done" = dashboardLabel("Done", status = "success")
    )
}

