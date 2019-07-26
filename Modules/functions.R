
#analysisAvail=c("16s18sits"=1,"rnaseq"=2,"metagenomic"=3,"gbs"=4)
analysisAvail=c("16s18sits"=1,"metagenomic"=3)

pipelines<-list()
for(module in names(analysisAvail)){
  source(paste0("Classes/",module,"_object.R"))
  pipelines[[module]]<-getpipeline()
}
esclavoHome<-getwd()
projectName<-reactiveVal("")
analysisType<-reactiveVal("")
projectsFolder="/home/sandro/Programas/ESCLAVO/projects"
projectConf<-reactiveVal(data.frame())

makeProjectDescription<-function(projectFolder,fastqFolder,analysisType, aVersion, pStatus,pPercent,lStep){
  rnames<-c("pfolder","ffolder","fqpattern","analysis","aversion","created","status","pPercent","lastStep")
  nfastq<-length(list.files(fastqFolder,"*.fastq.gz"))
  fastqPattern<-ifelse(nfastq == 0,".fastq",".fastq.gz")

  df<-data.frame(value=c(projectFolder,fastqFolder,fastqPattern,
                                       analysisType,aVersion,
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

parseTimes<-function(timesVector,type){
  if(length(timesVector)==0){return(0)}
  times<-lapply(timesVector, function(x){
    times<-as.numeric(strsplit(x,":")[[1]])
    
    switch(type,
      "seconds" = {times[1]*60*60 + times[2]*60 + times[3]},
      "minutes" = {times[1]*60 + times[2] + times[3]/60},
      "hours" = {times[1] + times[2]/60 + times[3]/60/60}
    )
  })
  times<-round(sum(unlist(times)),2)
  paste(times,type)
}
