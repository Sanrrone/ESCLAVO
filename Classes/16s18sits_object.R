source("Classes/16s18sits_tabcontent.R")
source("Classes/16s18sits_module.R")

getpipeline<-function(){
  new("pipeline",
      id = "16s18sits",
      name = "16S-18S/ITS",
      version = 1.0,
      image = "images/16srrna.png",
      steps = list(
        statusb = list(folder="0-raw",
                       software="MultiQC",
                       version="0.0.0",
                       spath="/home/sandro/.local/bin/multiqc",
                       stepName="Reads before QC",
                       stepID="rbqc",
                       icon = icon("stream"),
                       iconHTML='fa fa-stream',
                       tabcontentUI=statusbUIm("statusbmodule","step_rbqc","0-raw","MultiQC","0.0.0","/home/sandro/.local/bin/multiqc"),
                       tabcontentSrv=list(server=statusbTabModule,
                                          id="statusbmodule")),
        qc = list(folder="1-qc",
                  software="Prinseq",
                  version="0.0.0",
                  stepName="QC status",
                  stepID="qc",
                  icon = icon("filter"),
                  iconHTML='fa fa-filter',
                  tabcontentUI=qcTab("step_qc"),
                  tabcontentSrv=list(server=statusbTabModule,
                                     id="statusbmodule")),
        statusa = list(folder="1-qc",
                       software="Prinseq",
                       version="0.0.0",
                       stepName="Reads after QC",
                       stepID="raqc",
                       icon = icon("stream"),
                       iconHTML='fa fa-stream',
                       tabcontentUI=statusaTab("step_raqc"),
                       tabcontentSrv=list(server=statusbTabModule,
                                          id="statusbmodule")),
        taxInsight = list(folder="2-taxInsight",
                          software="dada2",
                          version="0.0.0",
                          stepName="Taxonomic counts",
                          stepID="tc",
                          icon = icon("microscope"),
                          iconHTML='fa fa-microscope',
                          tabcontentUI=taxCountTab("step_tc"),
                          tabcontentSrv=list(server=statusbTabModule,
                                             id="statusbmodule"))
      ),
      stepsOptions = list(
        color = "danger",
        collapsed = T
      )
  )
}