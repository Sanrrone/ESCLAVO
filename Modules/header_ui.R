# Module UI function
headerUIm <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  dashboardHeaderPlus(
    title = tagList(
      span(class = "logo-lg", "ESCLAVO"), 
      icon("heart")),
    enable_rightsidebar = T,
    rightSidebarIcon = "gears",
    left_menu = tagList(
      uiOutput("header_dropdown_proyectName"),
      uiOutput("header_dropdown_proyectType")
    )
  )
}