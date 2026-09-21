###############################################################
#  SHINY DASHBOARD 
###############################################################

install.packages(c("shinydashboard","DT","bslib"))


library(shiny)
library(shinydashboard)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(TH.data)
library(DT)

###############################
# Données
###############################
data("GBSG2", package="TH.data")

df <- GBSG2 %>%
  mutate(
    horTh = factor(horTh),
    menostat = factor(menostat),
    tgrade = factor(tgrade, ordered = TRUE),
    cens = factor(cens, labels = c("Censuré","Événement"))
  )

num_vars <- c("age","tsize","pnodes","progrec","estrec","time")
cat_vars <- c("horTh","menostat","tgrade","cens")

###############################
# UI
###############################
ui <- dashboardPage(
  
  dashboardHeader(title = "Analyse de survie"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Accueil", tabName = "overview", icon = icon("chart-pie")),
      menuItem("Analyse descriptive", tabName = "desc", icon = icon("table")),
      menuItem("Analyse de survie", tabName = "survival", icon = icon("heartbeat")),
      menuItem("Modèle de Cox", tabName = "cox", icon = icon("chart-line"))
    )
  ),
  
  dashboardBody(
    
    tabItems(
      
      ###############################################################
      # VUE D'ENSEMBLE
      ###############################################################
      tabItem(
        tabName = "overview",
        
        fluidRow(
          valueBoxOutput("n_patients"),
          valueBoxOutput("n_events"),
          valueBoxOutput("median_surv")
        ),
        
        box(
          title = "Aperçu du jeu de données",
          width = 12,
          status = "primary",
          solidHeader = TRUE,
          DTOutput("head_table")
        )
      ),
      
      ###############################################################
      # ANALYSE DESCRIPTIVE
      ###############################################################
      tabItem(
        tabName = "desc",
        
        fluidRow(
          box(
            title = "Filtres",
            width = 4,
            status = "info",
            solidHeader = TRUE,
            
            checkboxGroupInput(
              "filter_horTh",
              "Traitement hormonal",
              choices = levels(df$horTh),
              selected = levels(df$horTh)
            ),
            
            checkboxGroupInput(
              "filter_meno",
              "Statut ménopausique",
              choices = levels(df$menostat),
              selected = levels(df$menostat)
            ),
            
            selectInput("num_var","Variable numérique",choices=num_vars),
            selectInput("cat_var","Variable qualitative",choices=cat_vars)
          ),
          
          box(
            title = "Statistiques descriptives",
            width = 8,
            status = "primary",
            solidHeader = TRUE,
            verbatimTextOutput("desc_stats")
          )
        ),
        
        fluidRow(
          box(title="Histogramme", width=6, plotOutput("hist_plot")),
          box(title="Boxplot", width=6, plotOutput("box_plot"))
        ),
        
        fluidRow(
          box(title="Variable qualitative", width=12, plotOutput("bar_plot"))
        )
      ),
      
      ###############################################################
      # SURVIE
      ###############################################################
      tabItem(
        tabName = "survival",
        
        box(
          title = "Courbes de survie Kaplan–Meier",
          width = 12,
          status = "danger",
          solidHeader = TRUE,
          plotOutput("km_plot")
        )
      ),
      
      ###############################################################
      # COX
      ###############################################################
      tabItem(
        tabName = "cox",
        
        fluidRow(
          box(
            title = "Résumé du modèle de Cox",
            width = 6,
            status = "warning",
            solidHeader = TRUE,
            verbatimTextOutput("cox_summary")
          ),
          
          box(
            title = "Forest plot",
            width = 6,
            status = "warning",
            solidHeader = TRUE,
            plotOutput("cox_forest")
          )
        )
      )
    )
  )
)

###############################
# SERVER
###############################
server <- function(input, output) {
  
  df_filtered <- reactive({
    df %>%
      filter(
        horTh %in% input$filter_horTh,
        menostat %in% input$filter_meno
      )
  })
  
  ###############################
  # KPI
  ###############################
  
  output$n_patients <- renderValueBox({
    valueBox(nrow(df_filtered()), "Patientes", icon = icon("users"), color = "blue")
  })
  
  output$n_events <- renderValueBox({
    valueBox(sum(df_filtered()$cens=="Événement"),
             "Événements", icon = icon("heartbeat"), color = "red")
  })
  
  output$median_surv <- renderValueBox({
    fit <- survfit(Surv(time, cens=="Événement") ~ 1, data=df_filtered())
    valueBox(round(summary(fit)$table["median"],1),
             "Survie médiane", icon = icon("clock"), color = "green")
  })
  
  ###############################
  # TABLE
  ###############################
  
  output$head_table <- renderDT({
    datatable(head(df_filtered(),10), options = list(pageLength=5))
  })
  
  ###############################
  # DESCRIPTIF
  ###############################
  
  output$desc_stats <- renderPrint({
    summary(df_filtered()[,num_vars])
  })
  
  output$hist_plot <- renderPlot({
    ggplot(df_filtered(), aes_string(input$num_var)) +
      geom_histogram(fill="#1F78B4", bins=30, alpha=0.8) +
      theme_minimal()
  })
  
  output$box_plot <- renderPlot({
    ggplot(df_filtered(),
           aes_string(x="horTh", y=input$num_var, fill="horTh")) +
      geom_boxplot(alpha=0.8) +
      theme_minimal()
  })
  
  output$bar_plot <- renderPlot({
    ggplot(df_filtered(), aes_string(input$cat_var, fill=input$cat_var)) +
      geom_bar(alpha=0.85) +
      theme_minimal()
  })
  
  ###############################
  # SURVIE
  ###############################
  
  output$km_plot <- renderPlot({
    fit <- survfit(Surv(time, cens=="Événement") ~ horTh, data=df_filtered())
    
    ggsurvplot(
      fit,
      data=df_filtered(),
      pval=TRUE,
      risk.table=TRUE,
      conf.int=TRUE,
      ggtheme=theme_minimal(),
      palette=c("#E64B35","#4DBBD5")
    )$plot
  })
  
  ###############################
  # COX
  ###############################
  
  output$cox_summary <- renderPrint({
    model <- coxph(Surv(time, cens=="Événement") ~ age + horTh + tgrade + menostat,
                   data=df_filtered())
    summary(model)
  })
  
  output$cox_forest <- renderPlot({
    model <- coxph(Surv(time, cens=="Événement") ~ age + horTh + tgrade + menostat,
                   data=df_filtered())
    ggforest(model, data=df_filtered())
  })
}

shinyApp(ui, server)

