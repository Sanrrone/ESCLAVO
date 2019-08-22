if('rmarkdown' %in% rownames(installed.packages()) == FALSE) {install.packages('rmarkdown')}
if('rmdformats' %in% rownames(installed.packages()) == FALSE) {install.packages('rmdformats')}
if('DT' %in% rownames(installed.packages()) == FALSE) {install.packages('DT')}
if('dplyr' %in% rownames(installed.packages()) == FALSE) {install.packages('dplyr')}
if('reshape2' %in% rownames(installed.packages()) == FALSE) {install.packages('reshape2')}
if('plotly' %in% rownames(installed.packages()) == FALSE) {install.packages('plotly')}
if('stringr' %in% rownames(installed.packages()) == FALSE) {install.packages('stringr')}

library(rmarkdown)
library(rmdformats)
library(DT)
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(stringr)
options(scipen=999)

reportParams<-list()
############################## 0-raw
reportParams[["nfastq"]]<-length(list.files("*.fastq.gz",path = "FASTQFOLDER"))/2
allfiles<-list.files("PFOLDER",recursive =T)
#to remove report R files
reportParams[["allfiles"]]<-allfiles[unlist(lapply(allfiles,function(x)!grepl("^report.R[md]",x)))] 
reportParams[["proyectname"]]<-"PNAME"

############################## 1-qc
reportParams[["qcreport"]]<-"MULTIQCHTML"
qctable<-read.table("PFOLDER/1-qc/qc_summary.tsv",sep = "\t",header = T)
qctable$Sample<-rownames(qctable)
rownames(qctable)<-1:nrow(qctable)
qctable<-qctable[,c(ncol(qctable),1:(ncol(qctable)-1))]
reportParams[["qctable"]]<-qctable
#execution pdf

############################## 2-taxInsight
SummaryAbu<-read.table("ABUNDANCEFILE",sep = "\t",header = T,stringsAsFactors = F) #abundance file to include in report
tmp<-data.frame()
for(col in colnames(SummaryAbu)[-1:-7]){
  top10<-SummaryAbu[order(SummaryAbu[,col],decreasing = T)[1:10],]
  tmp<-merge(tmp,top10,all=T)
}
SummaryAbu<-tmp[,colnames(tmp)[-1:-4]]

reportParams[["SummaryAbu"]]<-SummaryAbu
reportParams[["taxcountFolder"]]<-"TAXFOLDER" #tax folder where pcoa.png is located.




reportParams[["outtype"]]<-"pdf"
render(input = "report.Rmd", params = reportParams, output_format = "pdf_document")

#executuon html
reportParams[["outtype"]]<-"html"
render(input = "report.Rmd", params = reportParams)


