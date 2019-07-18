source("Classes/16s18sits_tabcontent.R")
getpipeline<-function(){
  new("pipeline",
      id = "16s18sits",
      name = "16S-18S/ITS",
      version = 1.0,
      image = "images/16srrna.png",
      steps = list(
        statusb = list(folder="0-raw",
                       software="FastQC",
                       version="0.0.0",
                       stepName="Reads before QC",
                       stepID="rbqc",
                       icon = icon("stream"),
                       iconHTML='fa fa-stream',
                       tabcontent=statusbTab("step_rbqc")),
        qc = list(folder="1-qc",
                  software="Prinseq",
                  version="0.0.0",
                  stepName="QC status",
                  stepID="qc",
                  icon = icon("filter"),
                  iconHTML='fa fa-filter',
                  tabcontent=qcTab("step_qc")),
        statusa = list(folder="1-qc",
                       software="Prinseq",
                       version="0.0.0",
                       stepName="Reads after QC",
                       stepID="raqc",
                       icon = icon("stream"),
                       iconHTML='fa fa-stream',
                       tabcontent=statusaTab("step_raqc")),
        taxInsight = list(folder="2-taxInsight",
                          software="dada2",
                          version="0.0.0",
                          stepName="Taxonomic counts",
                          stepID="tc",
                          icon = icon("microscope"),
                          iconHTML='fa fa-microscope',
                          tabcontent=taxCountTab("step_tc"))
      ),
      stepsOptions = list(
        color = "danger",
        collapsed = T
      )
  )
}