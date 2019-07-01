# Module UI function
headerUIm <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  dashboardHeaderPlus(
    title = "ESCLAVO v1.0",
    enable_rightsidebar = F,
    rightSidebarIcon = "gears",
    left_menu = tagList(
      uiOutput("header_dropdown_proyectName"),
      uiOutput("header_dropdown_proyectType")
    )
  )
}