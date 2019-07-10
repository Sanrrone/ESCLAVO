
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

makeProjectDescription<-function(projectFolder, fastqFolder, analysisType, aVersion, pStatus){
  df<-data.frame(pfolder=projectFolder,
                 ffolder=fastqFolder,
                 analysis=analysisType,
                 aversion=aVersion,
                 created=Sys.Date(),
                 status=pStatus,
                 
                 stringsAsFactors = F)
  
  
  return(df)
}