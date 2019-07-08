getpipeline<-function(){
new("pipeline",
    id = "metagenomic",
    name = "Metagenomic",
    version = 1.0,
    image = "images/metagenomic.jpg",
    steps = list(
      statusb = "fastqc",
      qc = "prinseq",
      statusa = "fastqc",
      align = "centrifuge"
    ),
    stepsOptions = list(
      color = "danger",
      collapsed = T
    )
)
}