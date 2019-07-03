rsidebarModule<-function(input, output, session) {

  output$rsidebarcontent <- renderUI({
    content<-list()
    for(module in names(analysisAvail)){
      switch(module,
             "16s18sits" =  {
               content[["16s18sits"]]<-rightSidebarTabContent(
                 id = "16s18sits",
                 title = "16S-18S/ITS",
                 icon = "vial"
               )
             },
             "rnaseq" = {
               content[["rnaseq"]]<-rightSidebarTabContent(
                 id = "rnaseq",
                 title = "RNA-Seq",
                 icon="vials"
               )
             },
             "metagenomic" = {
               content[["metagenomic"]]<-rightSidebarTabContent(
                 id = "metagenomic",
                 icon = "dna",
                 title = "Metagenomic",
                 numericInput("obs", "Observations:", 10, min = 1, max = 100)
               )
             },
             "gbs" = {
               content[["gbs"]]<-rightSidebarTabContent(
                 id = "gbs",
                 icon = "sliders",
                 title = "Genotype by sequencing"
               )
             }
             
      )
    }
    
    rightSidebar( background = "dark", content)

  })
}