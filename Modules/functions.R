#analysisAvail=c("16s18sits"=1,"rnaseq"=2,"metagenomic"=3,"gbs"=4)
analysisAvail=c("16s18sits"=1,"metagenomic"=3)

pipelines<-list()
for(module in names(analysisAvail)){
  source(paste0("Classes/",module,"_object.R"))
  pipelines[[module]]<-getpipeline()
}


