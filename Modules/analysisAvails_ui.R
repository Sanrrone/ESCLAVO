# Module UI function
analysisAvailUIm <- function(id,analysisAvail) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  uiOutput(ns("analysisImages"))
  
}