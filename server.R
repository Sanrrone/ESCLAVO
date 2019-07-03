source("Modules/rsidebar_server.R")
function(input, output) { 

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
                 icon = "align-center",
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
    
    return(rightSidebar( background = "dark", content))
    
  })
output$itemListUI <- renderUI({
  userlist<-list()
  for(module in names(analysisAvail)){
    switch(module,
           "16s18sits" =  {
             userlist[["16s18sits"]]<-userListItem(
               src = dataURI(file="images/16srrna.png", mime="image/png"), 
               user_name = "16S, 18S/ITS", 
               description = "Updated: August 2019"
             )
           },
           "rnaseq" = {
             userlist[["rnaseq"]]<-userListItem(
               src = dataURI(file="images/rnaseq.jpg", mime="image/png"), 
               user_name = "RNA-Seq", 
               description = "Updated: August 2019"
             )
           },
           "metagenomic" = {
             userlist[["metagenomic"]]<-userListItem(
               src = dataURI(file="images/metagenomic.jpg", mime="image/png"), 
               user_name = "Metagenomic", 
               description = "Updated: August 2019"
             )
           },
           "gbs" = {
             userlist[["gbs"]]<-userListItem(
               src = dataURI(file="images/gbs.jpg", mime="image/png"), 
               user_name = "Genotype by sequencing", 
               description = "Updated: August 2019"
             )
           }
    )
  }
    
   return(userList(userlist))

})


  
}
