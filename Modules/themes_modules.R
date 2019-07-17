# Ui functions ------------------------------------------------------------
uiChangeThemeDropdown <- function(dropDownLabel = "Dashboard Theme")
{
  changeThemeChoices <- c(
    "Default" = NA,
    "Blue gradient" = "blue_gradient",
    "BoE website" = "boe_website",
    "Grey light" = "grey_light",
    "Grey dark" = "grey_dark",
    "OneNote" = "onenote",
    "Poor man's Flatly" = "poor_mans_flatly",
    "Purple gradient" = "purple_gradient"
  )
  
  ns <- NS("moduleChangeTheme")
  dropdownBlock(
    id = "mydropdown",
    title = "Dashboard Theme",
    icon = "sliders",
    selectizeInput(
      inputId = ns("dbxChangeTheme"),
      label = dropDownLabel,
      choices = changeThemeChoices
    )
  )

}

uiChangeThemeOutput <- function()
{
  ns <- NS("moduleChangeTheme")
  themeOutput <- tagList(
    uiOutput(ns("uiChangeTheme"))
  )
  
  return(themeOutput)
}


# Server functions --------------------------------------------------------
serverChangeTheme <- function(input, output, session)
{
  observeEvent(input$dbxChangeTheme,{
      output$uiChangeTheme <- renderUI({
        shinyDashboardThemes(theme = input$dbxChangeTheme)
      })
    }
  )
}