# Module UI function
rightsidebarUIm <- function(id,analysisAvail) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  rightSidebar( background = "dark",
                rightSidebarTabContent(
                  id = "esclavoOptions",
                  title = "Analysis options",
                  icon = "vials",
                  active = T,
                  uiOutput(ns("rsidebarcontent"))
                )
                
                
  )
  
}