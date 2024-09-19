#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(readxl)
library(ggplot2)
library(plyr)
library(shinythemes)
library(shinysky)
library(shinyBS)
library(dplyr)

# Define UI for application that takes collaborative cross 
ui <- fluidPage(theme=shinytheme("cerulean"),
  titlePanel("Comparative Collaborative Cross Proteomics Expression"),
  sidebarLayout(
    sidebarPanel(
      h2("Setup",
         bsButton("q1", label = "", icon("question"), style = "info", size = "extra-small"),
         bsPopover(id="q1", title = "Conlon Lab Collaborative Cross Data"),
         content = "This is not a real button",
         placement = "right", 
         trigger = "click",
         options = list(container = "body")
         ),
      checkboxGroupInput(inputId = "strains",
                         label = "Strains to Include:", 
                         choices = c("PWK", "C57", "AJ", "129", "WSB", "CAS", "NZO", "NOD"),
                         selected = c("PWK", "C57", "AJ", "129", "WSB", "CAS", "NZO", "NOD"),
                         inline = FALSE,
                         width = NULL,
                         choiceNames = NULL,
                         choiceValues = NULL
      ),
      selectInput(inputId = "normstrain", label = "Reference Strain", choices = c("C57", "PWK", "AJ", "129", "WSB", "CAS", "NZO", "NOD")),
      textInput(inputId = "proteinoi", 
                label = "Protein of Interest:", 
                value = "", width = NULL, placeholder = NULL
      ),
      shiny::actionButton(inputId = "go", class = "btn btn-lg btn-primary", label = "Setup Complete")
    ),
    mainPanel(
      tabsetPanel(id="alltabs", type = "tabs",
        tabPanel("Protein Abundance",
                 shiny::uiOutput("Setup"),
                 plotOutput("tab1expressionboxplot"),
                 downloadButton("DownloadPlot", "Download Plot")
                 ),
        tabPanel("Male/Female Ratio",
                 plotOutput("tab2expressionbarplot"),
                 downloadButton("DownloadPlot2", "Download Plot"),
                 plotOutput("tab3plot"),
                 downloadButton("DownloadPlot3", "Download Plot"),
                 tableOutput("tab3table"),
                 downloadButton("DownloadTable", "Download Table")
                 ),
                ),
    )
  )
)

server <- function(input, output, session) {

  df=NULL
  df2=NULL
  colors=c("#F0E442", "#555555", "#009E73", "#0072B2", "#56B4E9", "#D55E00", "#CC79A7", "#E69F00" )
  names(colors)<- c("AJ", "C57", "CAS", "NOD", "NZO", "PWK", "WSB", "129" )
  
  
  get_metadata<- reactive({
    proteinstouse<- input$proteinoi
    proteinstouse<-unlist(strsplit(proteinstouse, ","))
    return(proteinstouse)
  })
  get_strains<- reactive({
    strainstouse<- input$strains
    strainstouse<-unlist(strsplit(strainstouse, ","))
    return(strainstouse)
  })
  get_normalized_strain<- reactive({
    norm<- input$normstrain
    return(norm)
  })
  
  #Load in the collaborative cross data frame
  cc_df <- read.table("https://github.com/drlaurenwasson/ConlonLab-CollabX/raw/main/files/cc_df.txt", sep = "\t", row.names = 1)
  genes<- rownames(cc_df)
  cc_df2<- read.table("https://github.com/drlaurenwasson/ConlonLab-CollabX/raw/main/files/cc_df2.txt", sep = "\t")
  colnames(cc_df2)[4]<- "129"
  rownames(cc_df2)<- genes
  cc_df3<- read.table("https://github.com/drlaurenwasson/ConlonLab-CollabX/raw/main/files/cc_df3.txt", sep = "\t")
  colnames(cc_df3)[4]<- "129"
  rownames(cc_df3) = genes
  #Get a data frame for the proteins of interest for the strains of interest
  getdata<- function(protein,strain){
    allstrains<- c("PWK", "C57", "AJ", "129", "WSB", "CAS", "NZO", "NOD")
    protein<- protein[protein %in% rownames(cc_df)]
    df<- as.data.frame(t(cc_df[rownames(cc_df)==protein,]))
    df$strain<- c(rep("PWK",4), rep("C57",4), rep("AJ",4), rep("129",4), rep("WSB",4), rep("CAS",4), rep("NZO",4), rep("NOD",4))
    df$protein<- protein
    colnames(df)[1]<-"Expression"
    #Get all of the average expressions for each strain
    avg<- c(mean(df$Expression[1:4]), mean(df$Expression[5:8]), mean(df$Expression[9:12]), mean(df$Expression[13:16]), mean(df$Expression[17:20]), mean(df$Expression[21:24]), mean(df$Expression[25:28]), mean(df$Expression[29:32]))
    names(avg)<- unique(df$strain)
    df$sex<- c("M", "F", "M", "F", "F", "F", "M",
               "M", "F", "M", "F", "M", "M", "M",
               "F", "F", "F", "M", "F", "M", "F", 
               "M", "F", "M", "M", "F", "M", "F",
               "M", "F", "M", "F")
    df$ss<- paste(df$strain,df$sex)
    return(df)
  }
  
  
  #When "Setup Complete" is pressed:
  shiny::observeEvent(input$go,{
    #Update the tabs list
    showTab(inputId = "alltabs", target = "Protein Expression")
    #Update the proteins list
    proteinsupdate<- get_metadata()
    output$proteinstouse<- shiny::renderText({proteinsupdate})
    #Update the strains to include
    strainsupdate<- get_strains()
    output$strainstouse<- shiny::renderText({strainsupdate})
    #Update the strain to normalize to
    normalizeupdate<- get_normalized_strain()
    output$normalizedstrain<- shiny::renderText({normalizeupdate})
    #Get the data for the genes we want and the strains we want
    
    for (p in {proteinsupdate}){
      #Check to see if protein is in the database
      if (p %in% rownames(cc_df)){
        a<- getdata(p,{strainsupdate})
        df=rbind(df,a)
      }
      else{
        output$proteinstouse<- shiny::renderText("")
      }
      output$dataframe<- shiny::renderTable({
        df}, rownames = TRUE
      )
    }
    #Tab 1: Plot protein expression values
    #Subset only the strains you want
    df<- df[df$strain %in% {strainsupdate},]
    #Plot the box plot
    p<- ggplot(df, aes(x=strain, y=Expression, group = ss, fill = ss)) + 
      geom_boxplot() + labs(title=paste0("Normalized Protein Expression- ", sort(unique(df$protein))), x="Strain", y = "Avg. Abundance (Normalized)")+theme_classic() 
    output$tab1expressionboxplot <- renderPlot({
      plot(p)
    })
    
    output$DownloadPlot = downloadHandler(
      file= paste({proteinsupdate}, ".png", sep=""),
      content = function(file){
        ggsave(p, filename = file)
      }
    )
    
    #Tab 2: Plot male to female ratio
    df2<- as.data.frame(t(cc_df2[rownames(cc_df2)=={proteinsupdate},]))
    df2$strain<- colnames(cc_df2)
    df2$strain[4]<- "129"
    df2$strain<- factor(df2$strain, levels=c("AJ","C57","129","NOD", "NZO", "CAS", "PWK", "WSB"))
    colnames(df2)[1]<-"Log2FC"
    df2$color<- c("#D55E00", "#555555", "#F0E442", "#E69F00", "#CC79A7", "#009E73", "#56B4E9", "#0072B2" )
    #Subset the rows you want
    df2<- df2[df2$strain %in% {strainsupdate},]
    p2<- ggplot(df2,aes(x=strain, y=Log2FC))+
      geom_bar(stat="identity", fill = df2$color ) 
    p2<- p2+labs(title=paste0("Normalized Protein Expression- ", {proteinsupdate}), x="Strain", y = "Log2FC M/F")+
      theme_classic() 
    output$tab2expressionbarplot <- renderPlot({
      plot(p2)
    })
    
    output$DownloadPlot2 = downloadHandler(
      file= paste({proteinsupdate}, "log2.png", sep=""),
      content = function(file){
        ggsave(p2, filename = file)
      }
    )
    
    #Tab 3: Normalize to the reference strain2
    df3<- as.data.frame(t(cc_df3[rownames(cc_df3)=={proteinsupdate},]))
    df3$strain<- colnames(cc_df3)
    df3$strain[4]<- "129"
    df3$strain<- factor(df3$strain, levels=c("AJ","C57","129","NOD", "NZO", "CAS", "PWK", "WSB"))
    colnames(df3)[1]<-"Log2FC"
    df3$color<- c("#D55E00", "#555555", "#F0E442", "#E69F00", "#CC79A7", "#009E73", "#56B4E9", "#0072B2" )
    df3$adj<- df3$Log2FC/df3[rownames(df3)[rownames(df3) %in% {normalizeupdate}],1]
    #Subset the rows you want
    df3<- df3[df3$strain %in% {strainsupdate},]
    p3<- ggplot(df3,aes(x=strain, y=adj))+
      geom_bar(stat="identity", fill = df3$color ) 
    p3<- p3+labs(title=paste0("Adjusted Normalized Protein Expression- ", {proteinsupdate}), x="Strain", y = "Adj Log2FC M/F")+
      theme_classic() 
    output$tab3plot<- renderPlot({
      plot(p3)
    })
    
    output$DownloadPlot3 = downloadHandler(
      file= paste({proteinsupdate}, "adjusted.png", sep=""),
      content = function(file){
        ggsave(p3, filename = file)
      }
    )
    output$tab3table<- renderTable({df3})
    
    output$DownloadTable = downloadHandler(
      file= paste({proteinsupdate}, ".txt", sep=""),
      content = function(file){
        write.table(df3, filename = file)
      }
    )
  })

}

# Run the application 
shinyApp(ui = ui, server = server)
