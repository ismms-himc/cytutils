library(shiny)
library(shinydashboard)
library(shinyFiles)

fcs_channel_rename_tab_item <- tabItem(tabName = "fcs_channel_rename",
              h2("FCS Channel Renaming",
                style = "padding-left: 15px"),
              sidebarPanel(
                  p(paste(
                    "Please select the directory containing the FCS files requiring",
                    "renaming of channel names and/or descriptions.")),
                  shinyDirButton(id="fcs_channel_rename_dir", label="Choose directory", title="Select directory" )
                ),
                mainPanel(
                  htmlOutput("new_channel_rename_csv_generation_status"),
                  htmlOutput("import_channel_rename_csv_and_export_fcs_files_status")
                )
           )