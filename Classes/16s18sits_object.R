getpipeline<-function(){
  new("pipeline",
      id = "16s18sits",
      name = "16S-18S/ITS",
      image = "images/16srrna.png",
      steps = list(
        statusb = "fastqc",
        qc = "prinseq",
        merge = "pear",
        statusa = "fastqc",
        taxInsight = "dada2"
      ),
      stepsOptions = list(
        color = "danger",
        collapsed = T
      )
  )
}