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
      selectInput(inputId = "normstrain", label = "Reference Strain", choices = c("PWK", "C57", "AJ", "129", "WSB", "CAS", "NZO", "NOD")),
      textInput(inputId = "proteinoi", 
                label = "Protein(s) of Interest: (comma separated)", 
                value = "", width = NULL, placeholder = NULL
      ),
      shiny::actionButton(inputId = "go", class = "btn btn-lg btn-primary", label = "Setup Complete")
    ),
    mainPanel(
      #Create Tabs
      tabsetPanel(id="expression", type = "tabs",
                  tabPanel("Protein Expression Values"),
                  # Allow input information from the "Setup" panel if necessary.
                  shiny::uiOutput("Setup"),
                  plotOutput(outputId = "scatterplot"),
                  downloadButton("DownloadPlot", "Download Plot")
                  ),
    )
  )
)

server <- function(input, output, session) {

  get_metadata<- reactive({
    proteinstouse<- input$proteinoi
    proteinstouse<-unlist(strsplit(proteinstouse, ","))
    #print(proteinstouse)
    return(proteinstouse)
  })
  get_strains<- reactive({
    strainstouse<- input$strains
    strainstouse<-unlist(strsplit(strainstouse, ","))
    #print(strainstouse)
    return(strainstouse)
  })
  get_normalized_strain<- reactive({
    norm<- input$normstrain
    #print(norm)
    return(norm)
  })
  
  #df2=NULL
  #df4=NULL
  #dfrp=NULL
  
  #When "Setup Complete" is pressed:
  shiny::observeEvent(input$go,{
    #Update the proteins list
    proteinsupdate<- get_metadata()
    output$proteinstouse<- shiny::renderText({proteinsupdate})
    #Update the strains to include
    strainsupdate<- get_strains()
    output$strainstouse<- shiny::renderText({strainsMupdate})
    #Update the time point to normalize to
    normalizeupdate<- get_normalized_strain()
    output$normalizedstrain<- shiny::renderText({normalizeupdate})
    #Get the data for the genes we want and the strains we want
    
    for (p in {proteinsupdate}){
      #Check to see if protein is in the database
      if (p %in% rownames(rawvalues)){
        a<- getdata(p,{timepointsupdate},{normalizeupdate})
        df2=rbind(df2,a)
      }
      
    }
    #Subset only the timepoints you want
    df2<- df2[df2$day %in% {timepointsupdate},]
    #Plot the scatterplot
    output$scatterplot <- renderPlot({
      p<- ggplot(df2, aes(x=day, y=Normalized, group=protein, color=protein)) + 
        geom_line() +
        geom_point()+
        geom_errorbar(aes(ymin=Normalized-sd, ymax=Normalized+sd), width=.2,
                      position=position_dodge(0.05))
      p=p+labs(title="Protein Expression", x="Day", y = "Avg. Abundance (Normalized)")+
        theme_classic()
      plot(p)
      observeEvent(input$update, print(as.numeric(input$update)))
    })
    p<- ggplot(df2, aes(x=day, y=Normalized, group=protein, color=protein)) + 
      geom_line() +
      geom_point()+
      geom_errorbar(aes(ymin=Normalized-sd, ymax=Normalized+sd), width=.2,
                    position=position_dodge(0.05))
    p=p+labs(title="Protein Expression", x="Day", y = "Avg. Abundance (Normalized)")+
      theme_classic()
    
    output$DownloadPlot = downloadHandler(
      file= paste({proteinsupdate}, {timepointsupdate}, ".png", sep=""),
      content = function(file){
        ggsave(p, filename = file)
      }
    )
  })
  
  #Load in the collaborative cross data frame
  masterdf <- read.txt("https://github.com/drlaurenwasson/ConlonLab-MMExpression/raw/main/files/1-s2.0-S1534580723001818-mmc2.csv", row.names = 1)
  genes<- masterdf$Gene.Symbol
  rawvalues<- as.data.frame(masterdf[,29:52])
  rownames(rawvalues)<- make.names(genes, unique = TRUE)
  
  #Get a Protein data frame for the proteins of interest
  getdata<- function(protein,timepoint,normtp){
    tp<- c("E09.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5")
    protein<- protein[protein %in% rownames(rawvalues)]
    df<- as.data.frame(t(rawvalues[rownames(rawvalues)==protein,]))
    df$day<- c(rep("E09.5",3), rep("E10.5",3), rep("E11.5",3), rep("E12.5",3), rep("E13.5",3), rep("E14.5",3), rep("E15.5",3), rep("E16.5",3))
    df$protein<- protein
    colnames(df)[1]<-"Expression"
    #Get all of the average expressions for each time point
    avg<- c(mean(df$Expression[1:3]), mean(df$Expression[4:6]), mean(df$Expression[7:9]), mean(df$Expression[10:12]), mean(df$Expression[13:15]), mean(df$Expression[16:18]), mean(df$Expression[19:21]), mean(df$Expression[22:24]))
    names(avg)<- tp
    #Normalize to the time point we want
    df$Normalized<- df$Expression/(avg[names(avg) %in% normtp])
    
    data_summary <- function(data, varname, groupnames){
      require(plyr)
      summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
          sd = sd(x[[col]], na.rm=TRUE))
      }
      data_sum<-ddply(data, groupnames, .fun=summary_func,
                      varname)
      data_sum <- rename(data_sum, c("mean" = varname))
      return(data_sum)
    }
    df2 <- data_summary(df, varname="Normalized", 
                        groupnames=c("day", "protein"))
    head(df2)
    return(df2)
  }
  
  #RNA functions
  ## Get Metadata
  get_metadata_rna<- reactive({
    rnastouse<- input$rnaoi
    rnastouse<-unlist(strsplit(rnastouse, ","))
    #print(strainstouse)
    return(rnastouse)
  })
  get_timepoints_rna<- reactive({
    rnatimepointstouse <- input$rnatimepoint
    #print(timepointstouse)
    return(rnatimepointstouse)
  })
  get_normalized_timepoint_rna<- reactive({
    rnanormtimepoint<- input$normalize2
    #print(normtimepoint)
    return(rnanormtimepoint)
  })
  
  #Load in the RNA data frame
  rnavalues<- read.table("https://raw.githubusercontent.com/drlaurenwasson/ConlonLab-MMExpression/main/files/dds_combined_normalized_counts_proteingenes.txt", sep = "\t")
  #Get an RNA data frame for the RNAs of interest
  getdatarna<- function(gene, rtimepoint, rnormtp){
    rtp<- c("E09.5", "E10.5", "E12.5", "E14.5", "E16.5")
    gene<- gene[gene %in% rownames(rnavalues)]
    df3<- as.data.frame(t(rnavalues[rownames(rnavalues)==gene,]))
    df3$day<- c(rep("E09.5",3), rep("E10.5",3), rep("E12.5",3), rep("E14.5",3), rep("E16.5",3))
    df3$gene<- gene
    colnames(df3)[1]<-"Expression"
    #Get normalized data
    ravg<- c(mean(df3$Expression[1:3]), mean(df3$Expression[4:6]), mean(df3$Expression[7:9]), mean(df3$Expression[10:12]), mean(df3$Expression[13:15]))
    names(ravg)<- rtp
    #Normalize to the time point we want
    df3$Normalized<- df3$Expression/(ravg[names(ravg) %in% rnormtp])
    
    data_summary2 <- function(data, varname, groupnames){
      require(plyr)
      summary_func2 <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
          sd = sd(x[[col]], na.rm=TRUE))
      }
      data_sum2<-ddply(data, groupnames, .fun=summary_func2,
                       varname)
      data_sum2 <- rename(data_sum2, c("mean" = varname))
      return(data_sum2)
    }
    df4 <- data_summary2(df3, varname="Normalized", 
                         groupnames=c("day", "gene"))
    head(df4)
    return(df4)
  }
  
  #When "rnaupdate" is pressed
  shiny::observeEvent(input$rnaupdate,{
    #Update the RNA list
    rnaupdate<- get_metadata_rna()
    output$rnasstouse<- shiny::renderText({rnaupdate})
    #Update the time points to include
    rnatimepointsupdate<- get_timepoints_rna()
    output$rnatimepointstouse<- shiny::renderText({rnatimepointsupdate})
    #Update the time point to normalize to
    rnanormalizeupdate<- get_normalized_timepoint_rna()
    output$rnanormalizedtimepoint<- shiny::renderText({rnanormalizeupdate})
    
    #Get the data for the genes we want and the time points we want
    #A "for loop" iterates through all x in y. In this example, it will perform the code inside the {} for every gene in the variable 'genes' that I defined above.
    for (rna in {rnaupdate}){
      #Check to see if rna is in the database
      if (rna %in% rownames(rnavalues)){
        b<- getdatarna(rna,{rnatimepointsupdate},{rnanormalizeupdate})
        df4=rbind(df4,b)
      }
      
    }
    #Subset only the timepoints you want
    df4<- df4[df4$day %in% {rnatimepointsupdate},]
    #Plot the scatterplot
    output$rnascatterplot <- renderPlot({
      p2<- ggplot(df4, aes(x=day, y=Normalized, group=gene, color=gene)) + 
        geom_line() +
        geom_point() + 
        geom_errorbar(aes(ymin=Normalized-sd, ymax=Normalized+sd), width=.2,
                      position=position_dodge(0.05))
      p2<- p2+labs(title="RNA Expression", x="Day", y = "Avg. Abundance (Normalized)")+
        theme_classic() 
      plot(p2)
      observeEvent(input$rnaupdate, print(as.numeric(input$rnaupdate)))
    })
    
    
    output$rnaDownloadPlot = downloadHandler(
      file= paste({rnaupdate}, {rnatimepointsupdate}, ".png", sep=""),
      content = function(file){
        ggsave(p, filename = file)
      }
    )
  })
  
  #RNA + Protein Functions
  ## Get metadata
  get_metadata_rnaprotein<- reactive({
    rnaprotstouse<- input$rnaproteinoi
    rnaprotstouse<-unlist(strsplit(rnaprotstouse, ","))
    return(rnaprotstouse)
  })
  #Get the RNA time points
  get_timepoints_rnaproteinr<- reactive({
    rnaproteintimepointstouser <- input$rnaproteintimepointr
    return(rnaproteintimepointstouser)
  })
  #Get the protein time points
  get_timepoints_rnaproteinp<- reactive({
    rnaproteintimepointstousep <- input$rnaproteintimepointp
    return(rnaproteintimepointstousep)
  })
  #Get the normalized time point
  get_normalized_timepointrp<- reactive({
    normtimepoint<- input$normalize3
    return(normtimepoint)
  })
  #Get a data frame for the RNAs of interest
  getdatarnaprotein<- function(gene, rnatimepoint, proteintimepoint, normtpr){
    #Get the RNA
    rtp<- c("E09.5", "E10.5", "E12.5", "E14.5", "E16.5")
    gene<- gene[gene %in% rownames(rnavalues)]
    df3<- as.data.frame(t(rnavalues[rownames(rnavalues)==gene,]))
    df3$day<- c(rep("E09.5",3), rep("E10.5",3), rep("E12.5",3), rep("E14.5",3), rep("E16.5",3))
    df3$gene<- gene
    colnames(df3)[1]<-"Expression"
    #Get normalized data
    ravg<- c(mean(df3$Expression[1:3]), mean(df3$Expression[4:6]), mean(df3$Expression[7:9]), mean(df3$Expression[10:12]), mean(df3$Expression[13:15]))
    names(ravg)<- rtp
    #Normalize to the time point we want
    df3$Normalized<- df3$Expression/(ravg[names(ravg) %in% normtpr])
    data_summary2 <- function(data, varname, groupnames){
      require(plyr)
      summary_func2 <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
          sd = sd(x[[col]], na.rm=TRUE))
      }
      data_sum2<-ddply(data, groupnames, .fun=summary_func2,
                       varname)
      data_sum2 <- rename(data_sum2, c("mean" = varname))
      return(data_sum2)
    }
    df4 <- data_summary2(df3, varname="Normalized", 
                         groupnames=c("day", "gene"))
    #Get the Protein
    ptp<- c("E09.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5")
    normtpp<- normtpr
    protein=gene
    protein<- protein[protein %in% rownames(rawvalues)]
    df5<- as.data.frame(t(rawvalues[rownames(rawvalues)==protein,]))
    df5$day<- c(rep("E09.5",3), rep("E10.5",3), rep("E11.5",3), rep("E12.5",3), rep("E13.5",3), rep("E14.5",3), rep("E15.5",3), rep("E16.5",3))
    df5$gene<- protein
    colnames(df5)[1]<-"Expression"
    #Get all of the average expressions for each time point
    avg<- c(mean(df5$Expression[1:3]), mean(df5$Expression[4:6]), mean(df5$Expression[7:9]), mean(df5$Expression[10:12]), mean(df5$Expression[13:15]), mean(df5$Expression[16:18]), mean(df5$Expression[19:21]), mean(df5$Expression[22:24]))
    names(avg)<- ptp
    #Normalize to the time point we want
    df5$Normalized<- df5$Expression/(avg[names(avg) %in% normtpp])
    
    data_summary <- function(data, varname, groupnames){
      require(plyr)
      summary_func <- function(x, col){
        c(mean = mean(x[[col]], na.rm=TRUE),
          sd = sd(x[[col]], na.rm=TRUE))
      }
      data_sum<-ddply(data, groupnames, .fun=summary_func,
                      varname)
      data_sum <- plyr::rename(data_sum, c("mean" = varname))
      return(data_sum)
    }
    df6 <- data_summary(df5, varname="Normalized", 
                        groupnames=c("day", "gene"))
    df4$data<- "RNA"
    df6$data<- "Protein"
    df7<- rbind(df4,df6)
    return(df7)
  }
  
  
  #When "rnaproteinupdate" is pressed
  shiny::observeEvent(input$rnaprotupdate,{
    #Update the RNA list
    rnaprotupdater<- get_metadata_rnaprotein()
    output$rnasprotstouser<- shiny::renderText({rnaprotupdater})
    #Update the Protein time points to include
    rnaprottimepointsupdatep<- get_timepoints_rnaproteinp()
    output$rnaprottimepointstousep<- shiny::renderText({rnaprottimepointsupdatep})
    #Update the RNA time points to include
    rnaprottimepointsupdater<- get_timepoints_rnaproteinr()
    output$rnaprottimepointstouser<- shiny::renderText({rnaprottimepointsupdater})
    #Update the normalization time point
    rnaprotnormalizeupdate<- get_normalized_timepointrp()
    output$rnaprotnormalizedtimepoint<- shiny::renderText({rnaprotnormalizeupdate})
    
    #Get the data for the genes we want and the time points we want
    #A "for loop" iterates through all x in y. In this example, it will perform the code inside the {} for every gene in the variable 'genes' that I defined above.
    for (rna in {rnaprotupdater}){
      #Check to see if rna is in the database
      if (rna %in% rownames(rnavalues)){
        b<- getdatarnaprotein(rna,{rnaprottimepointsupdater}, {rnaprottimepointsupdatep},{rnaprotnormalizeupdate})
        dfrp=rbind(dfrp,b)
      }
      
    }
    #Subset only the timepoints you want
    dfrp$group<- paste0(dfrp$gene,dfrp$data)
    dfrp<- dfrp[dfrp$day %in% c({rnaprottimepointsupdater},{rnaprottimepointsupdatep}),]
    sf <- max(dfrp$Normalized)
    #Plot the scatterplot
    output$rnaprotscatter <- renderPlot({
      p3<- ggplot(dfrp) +
        geom_line(aes(x=day, y=Normalized, group=group, color=gene)) +
        geom_point(aes(x=day, y=Normalized, group=group, color=gene, shape = data),size = 2) + 
        scale_y_continuous(name = "RNA Expression", sec.axis = sec_axis( ~.*sf, name="Protein Expression")) + 
        guides(colour = guide_legend(override.aes = list(size=2)))
      p3<- p3+labs(title="RNA/Protein Expression")+
        theme_classic() 
      plot(p3)
      observeEvent(input$rnaprotupdate, print(as.numeric(input$rnaprotupdate)))
    })
    
    output$rnaprotDownloadPlot = downloadHandler(
      file= paste({rnaprotupdate}, {rnaprottimepointsupdatep}, {rnaprottimepointsupdater},".png", sep=""),
      content = function(file){
        ggsave(p, filename = file)
      }
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
