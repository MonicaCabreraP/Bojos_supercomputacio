#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(gtools)
library(stringr)

#setwd("/home/vruizser/PhD/2019-2020/Bojos_per_la_supercomputacio/shiny_app")
list_of_files = mixedsort(list.files(path= "./models",
                                     pattern = c("confusionMatrix_predictors_",".rds"), full.names = T))
cm_list = lapply(list_of_files, function(x) readRDS(x)) 
one= str_split(list_of_files, "predictors_", simplify = T)
one = as.data.frame(one)
splitted = str_split(str_remove(one$V2, ".rds"), "_")


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   title = "Machine learning - SKCM classifier",
   #titlePanel("Machine learning - SKCM classifier"),
   
   # Show a plot of the generated distribution
  plotOutput("distPlot", width = "100%", height = "600px"),
  
  hr(),
  
  # fluidRow(
  #   column(12, align = "center",
  #          h4("Most mutated genes in melanoma"),
  #          sliderInput(inputId = "bins",
  #                                         "Number of selected genes:",
  #                                         min = 5,
  #                                         max = 50,
  #                                         value = 10),
  #          div(style =  "font-size : 500px;")
  #         
  #   )
  checkboxGroupInput("options", "Variables to train the model:",
                     c("Treatment" = "Rx",
                       "Mutation" = "BRAF.V600.Mut",
                       "BMI" = "BMI",
                       "Stage" = "Stage",
                       "Sex" = "Sex", 
                       "LDH" = "LDH",
                       "Age" = "Age"),
                     selected = c("BRAF.V600.Mut"),
                     inline = T)
  )
   

      




# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$distPlot <- renderPlot({
      # # generate bins based on input$bins from ui.R
      # x    <- faithful[, 2] 
      # bins <- seq(min(x), max(x), length.out = input$bins + 1)
      # 
      # # draw the histogram with the specified number of bins
      # hist(x, breaks = bins, col = 'darkgray', border = 'white')
      input_sorted = input$options
      index = lapply(splitted,
                     function(x) identical (sort(x), input_sorted))
      cm = cm_list[unlist(index)][[1]]
      
      layout(matrix(c(1,1,2)))
      par(mar=c(2,2,2,2))
      plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
      title('CONFUSION MATRIX', cex.main=2)
      
      # create the matrix
      rect(150, 430, 240, 370, col='#3F97D0')
      text(195, 435, row.names(cm$table)[1], cex=2.5)
      rect(250, 430, 340, 370, col='#F7AD50')
      text(295, 435, row.names(cm$table)[2], cex=2.5)
      text(125, 370, 'Predicted', cex=2.5, srt=90, font=2)
      text(245, 450, 'Actual', cex=2.5, font=2)
      rect(150, 305, 240, 365, col='#F7AD50')
      rect(250, 305, 340, 365, col='#3F97D0')
      text(140, 400, row.names(cm$table)[1], cex=2.5, srt=90)
      text(140, 335, row.names(cm$table)[2], cex=2.5, srt=90)

      # add in the cm results
      res <- as.numeric(cm$table)
      text(195, 400, res[1], cex=4, font=2, col='white')
      text(195, 335, res[2], cex=4, font=2, col='white')
      text(295, 400, res[3], cex=4, font=2, col='white')
      text(295, 335, res[4], cex=4, font=2, col='white')
      # 
      # add in the specifics
      plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
      text(10, 85, names(cm$byClass[1]), cex=2.5, font=2)
      text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=2.5)
      text(30, 85, names(cm$byClass[2]), cex=2.5, font=2)
      text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=2.5)
      text(50, 85, names(cm$byClass[5]), cex=2.5, font=2)
      text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=2.5)
      text(70, 85, names(cm$byClass[6]), cex=2.5, font=2)
      text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=2.5)
      text(90, 85, names(cm$byClass[7]), cex=2.5, font=2)
      text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=2.5)

      # add in the accuracy information
      text(30, 35, names(cm$overall[1]), cex=3, font=2)
      text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=2.5)
      text(70, 35, names(cm$overall[2]), cex=3, font=2)
      text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=2.5)
      
      
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

