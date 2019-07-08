source("Modules/rsidebar_server.R")
source("Modules/analysisAvails_server.R")
source("Modules/mainButtons_server.R")

function(input, output) {

callModule(rsidebarModule,"rsidebarmodule")
callModule(analysisAvailsModule,"analysisavails")
callModule(mainButtonsModule,"mainButtons")


}
