
#analysisAvail=c("16s18sits"=1,"rnaseq"=2,"metagenomic"=3,"gbs"=4)
analysisAvail=c("16s18sits"=1)

pipelines<-list()
for(module in names(analysisAvail)){
  source(paste0("Classes/",module,"_object.R"))
  pipelines[[module]]<-getpipeline()
}
esclavoHome<-getwd()
projectName<-reactiveVal("")
analysisType<-reactiveVal("")
projectFolder<-reactiveVal("/home/sandro/Programas/ESCLAVO/projects")
projectConf<-reactiveVal(data.frame())

############################# set cpu usage value
cores<-as.numeric(system("nproc",wait = T,intern = T))
getcpudf<-function(){
  con <- file("/proc/stat","r");cpuvec1 <- readLines(con,n=cores+1);close(con)
  cpudf1<-strsplit(cpuvec1," |  ") %>% bind_cols()
  colnames(cpudf1)<-cpudf1[1,]
  cpudf1<-data.frame(cpudf1[-1,-1])
  cpudf1
}

#cpudf1<-reactiveVal(cpudf1)
##############################

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

getCPUusage<-function(cpudf1){
  cpudf2<-getcpudf()
  cpudf<-lapply(colnames(cpudf1), function(x){
    a<-as.numeric(cpudf1[,x])
    b<-as.numeric(cpudf2[,x])
    
    user<-b[1]
    nice<-b[2]
    system<-b[3]
    idle<-b[4]
    iowait<-b[5]
    irq<-b[6]
    softirq<-b[7]
    steal<-b[8]
    guest<-b[9]
    guest_nice<-b[10]
    
    prevuser<-a[1]
    prevnice<-a[2]
    prevsystem<-a[3]
    previdle<-a[4]
    previowait<-a[5]
    previrq<-a[6]
    prevsoftirq<-a[7]
    prevsteal<-a[8]
    prevguest<-a[9]
    prevguest_nice<-a[10]
    
    PrevIdle<-previdle + previowait
    Idle<-idle + iowait
    
    PrevNonIdle<-prevuser + prevnice + prevsystem + previrq + prevsoftirq + prevsteal
    NonIdle<-user + nice + system + irq + softirq + steal
    
    PrevTotal<-PrevIdle + PrevNonIdle
    Total<-Idle + NonIdle
    
    totald<-Total - PrevTotal
    idled<-Idle - PrevIdle
    
    usage<-(totald - idled)/totald*100
    

    #usage<- 100*((b[1]+b[2]+b[3]+b[5]+b[6]+b[7]) - (a[1]+a[2]+a[3]+a[5]+a[6]+a[7])) / ((b[1]+b[2]+b[3]+b[4]+b[5]+b[6]+b[7]) - (a[1]+a[2]+a[3]+a[4]+a[5]+a[6]+a[7]))
    usage<-ifelse(is.nan(usage),0,usage)
    
    data.frame(usage=round(usage,2), cpu=x, stringsAsFactors = F)
  }) %>% bind_rows()
  #cpudf1(cpudf2)
  cpudf
}

