library(shiny)
library(shinydashboard)
source("./cytof_toolkit_helper_functions.R")
source_dir("./ui_modules")

ui <- dashboardPage(
  dashboardHeader(title = "CyTOF QC App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Sample Background Tracking", 
                tabName = "sample_background_app",
                icon = icon("flask"))
    )
  ),
  dashboardBody(
    tabItems(
      sample_background_app_tab_item
    )
  )
)