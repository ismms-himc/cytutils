library(shiny)
library(shinydashboard)
source("./cytof_toolkit_helper_functions.R")
source_dir("./ui_modules")

ui <- dashboardPage(
  dashboardHeader(title = "CyTOF QC App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("CyTOF QC Report Generator", 
          tabName = "qc_report_app", 
          icon = icon("check-circle-o")),
      menuItem("Sample Background Tracking", 
                tabName = "sample_background_app",
                icon = icon("flask"))
    )
  ),
  dashboardBody(
    tabItems(
      qc_report_app_tab_item,
      sample_background_app_tab_item
    )
  )
)