library(shiny)
library(shinydashboard)
source("./helper_functions.R")
source_dir("./ui_modules")

ui <- dashboardPage(
  dashboardHeader(title = "Cytutils"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("FCS Channel Renaming", 
                tabName = "fcs_channel_rename", 
                icon = icon("wrench")),
      # generate2dJsDivergenceDataFrame
      menuItem("JS Divergence (Multi-Pair)", 
                tabName = "js_divergence_multi_pair",
                icon = icon("keyboard-o")),
      # calculate2dJsDivergence
      menuItem("JS Divergence (Single-Pair)", 
                tabName = "js_divergence_single_pair",
                icon = icon("calculator")),
      menuItem("AOF", 
          tabName = "aof",
          icon = icon("superpowers")),
      menuItem("Greedy AOF", 
          tabName = "greedy_aof",
          icon = icon("eercast"))
    )
  ),
  dashboardBody(
    # NB: The tabItems below are sourced from "./ui_modules"
    tabItems(
      fcs_channel_rename_tab_item
      # js_divergence_multi_pair_tab_item,
      # js_divergence_single_pair_tab_item,
      # aof_tab_item,
      # greedy_aof_tab_item
    )
  )
)