#Create shiny app which lets Users type in gene ID OR symbol and have validated the IDs
library(shiny)
library(bslib)
library(tidyverse)
library(rms)

z_score <- read.csv("data/z_scores_validated_id.csv", header = TRUE)
norm_1 <- read.csv("data/norm_counts_1.csv", header = TRUE)
norm_2 <- read.csv("data/norm_counts_2.csv", header = TRUE)
norm_3 <- read.csv("data/norm_counts_3.csv", header = TRUE)

norm_counts <- rbind(norm_1, norm_2, norm_3)

norm_counts <- norm_counts %>% arrange(sp)

ui <- page_sidebar(
  title = "Gene Trajectories",
  sidebar = sidebar(
    helpText(
      "IDs validated with FlyBase version FB2024_03"
    ),
    helpText(
      "If Gene Symbol has special characters, use FBgn Number."
    ),
    textInput("gene_id",
              label = "FBgn Number or Gene Symbol",
              value = "FBgn0024248"),
    selectInput(
      "tf",
      label="Choose Timeframe",
      choices = list(
        "Day",
        "Survival",
        "Sampling Point"
      ),
      selected = "Day"
    )
  ),
  card(textOutput("Id_status")),
  card(plotOutput("clus_plot")),
  card(plotOutput("gene_plot"))
)

server <- function(input, output){
  
  output$clus_plot <- renderPlot({
    if(input$gene_id %in% z_score$validated_id | input$gene_id %in% z_score$current_symbol){
        z_scores_clus <- z_score %>% filter(Official_ID == unique(z_score %>% filter(validated_id == input$gene_id | current_symbol == input$gene_id) %>% pull(Official_ID)))
     
      if(input$tf == "Day"){
        z_score_graph_st <- z_scores_clus %>% 
          ggplot(aes(x=day, y=Z_score))+xlab("Day")
      }else if (input$tf == "Survival"){
        z_score_graph_st <- z_scores_clus %>% 
          ggplot(aes(x=survivorship, y=Z_score))+scale_x_reverse()+xlab("Survival")
      }else{
        z_score_graph_st <- z_scores_clus %>% 
          ggplot(aes(x=sp, y=Z_score))+xlab("Sampling Point")
      }
      
      z_score_graph_st+
        geom_smooth(method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), se = FALSE, show.legend = FALSE,aes(group=gene_id),color="gray")+
        geom_smooth(method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), se = FALSE, show.legend = FALSE, color="black")+
        ggtitle(unique(z_scores_clus$Official_ID))+
        ylab("Z Scores")+
        theme_classic()
    }
  })
  
  output$gene_plot <- renderPlot({
    
    if(input$gene_id %in% norm_counts$validated_id | input$gene_id %in% norm_counts$current_symbol){
      gene_counts_long <- norm_counts %>% filter(validated_id == input$gene_id |current_symbol == input$gene_id)
      if(input$tf == "Day"){
        gene_counts_long$day <- as.numeric(gene_counts_long$day)
        gene_norm_graph <- gene_counts_long %>%
          ggplot(aes(x=day, y=Norm_Counts))+xlab("Day")
      }else if (input$tf == "Survival"){
        gene_norm_graph <- gene_counts_long %>%
          ggplot(aes(x=survivorship, y=Norm_Counts))+xlab("Survival")+scale_x_reverse()
      }else{
        gene_norm_graph <- gene_counts_long %>%
          ggplot(aes(x=sp,y=Norm_Counts))+xlab("Sampling Point")
      }
      
      title_id <- paste(unique(gene_counts_long$current_symbol)," (", unique(gene_counts_long$gene_id),")", sep = "")
      
      gene_norm_graph+
        geom_smooth(method = "lm", formula = y~rcs(x,quantile(x, c(0,0.25,0.5,0.75,1))), show.legend = FALSE, color = "red")+
        geom_point(size=2)+
        ylab("Normalized Read Counts")+
        ggtitle(title_id)+
        theme_classic()
    }
  })
  
  output$Id_status <- renderText({
    if(input$gene_id %in% norm_counts$validated_id | input$gene_id %in% norm_counts$current_symbol){
      norm_counts_report <- TRUE
    }else{
      norm_counts_report <- FALSE
    }
    if(input$gene_id %in% z_score$validated_id | input$gene_id %in% z_score$current_symbol){
      z_score_report <- TRUE
    }else{
      z_score_report <- FALSE
    }
    
    if(z_score_report & norm_counts_report){
      z_scores_gene <- z_score %>% filter(validated_id == input$gene_id |current_symbol == input$gene_id)
      paste(unique(z_scores_gene$current_symbol),"'s (", unique(z_scores_gene$validated_id),") expression chages with aging and is in ", unique(z_scores_gene$Official_ID),".",sep = "")
    }else if(z_score_report == FALSE & norm_counts_report){
      gene_counts_long <- norm_counts %>% filter(validated_id == input$gene_id |current_symbol == input$gene_id)
      paste(unique(gene_counts_long$current_symbol),"'s (", unique(gene_counts_long$validated_id),") expression does not significantly change with aging.", sep = "")
    }else if(Z_score_report == FALSE & norm_counts_report == FALSE){
      paste(input$gene_id, " was not identified in our analysis.", sep = "")
    }
  })
  
}

shinyApp(ui = ui, server =server)
