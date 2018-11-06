library(shiny)
library("V8")
library(shinyjs)
#library(MSnbase)
library(formattable)
library(stringr)
require(DT, quietly = TRUE) 
library(prozor)
library(markdown)

load("rule_db_multiple_sub_raw.RData")
source('helper.r')

textInputRow<-function (inputId, label, value = "") 
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}

shinyUI(navbarPage("MESSAR 0.1 (MEtabolite SubStructure Auto-Recommender)",
  
                   
    tabPanel("A) Start a run",
     shinyjs::useShinyjs(),
     shinyjs::extendShinyjs(text = "shinyjs.refresh = function() { location.reload(); }"),
             
     column(5,  
          
        br(),
        h4("Please paste your MS/MS spectrum into the field below:"), 
        textAreaInput("blank_file1", label = '',width=500,height=200),
        
        br(),
        
        h4("[Optional] Please paste the mass differences into the field below:"), 
        textAreaInput("blank_file2", label = '',width=500,height=150),
        
        textInput("prec_mz", h4("[Recommended] Precursor mass:"), value = "")),
  
    column(7,
        br(),
        
        numericInput("Relative", h4("Relative intensity threshold (base peak %)"),
                    min = 0, max = 99, value = 1, width = '500px'),
        br(),
        
        numericInput("ppm_search", h4("Tolerance [ppm] for masses and mass differences"),
                    min = 0, max = 50, value = 20, width = '500px'),
        br(),
        
        checkboxInput("fdr_control", label = "Filtering rules with a FDR cutoff at 0.05", value = TRUE, width = '500px'),
        br(),
        
        br(),
        tags$head(
          tags$style(HTML('#exampleButton{background-color:lightblue}'))
        ),
        actionButton("exampleButton", "Load example",style='padding:6px; font-size:150%'),
        br(),
        
        br(),
        tags$head(
          tags$style(HTML('#goButton{background-color:lightgreen}'))
        ),
        actionButton("goButton", "Submit",style='padding:6px; font-size:150%'),
        br(),
        
        br(),
        tags$head(
          tags$style(HTML('#killButton{background-color:orange}'))
        ),
        actionButton("killButton", "Clear",style='padding:6px; font-size:150%'),
        br(),
        
        br(),
        br(),
        em('Messages from the server:'),
        br(),
        br(),
        textOutput("blank_message1")
    )),

    tabPanel("B) Annotated features",
             
             tags$style("#blank_message2 {font-size:20px;
          color:red;
                        display:block; }"),        
             
       br(),
       div(style="display: inline-block;vertical-align:top; width: 550px;", uiOutput("blank_message2")),
       br(),     
       h4("Here is the list of annotated features (masses and mass differences):"), 
       br(),
       dataTableOutput("table1"),
       br(),
       downloadButton("annotated_rules", "Download matched rules",style='padding:6px; font-size:150%'),
       br(),
       plotOutput("plot_fdr",width = '1600px')
       ),

    tabPanel("C) Sub-structure suggestions",
       column(5,
       br(),
       selectInput("score_type", label= h4("Please select a scoring method:"), 
                   c("Lift sum [Most informative but less confident]"="L1", 
                     "Lift median [Most informative but less confident]"="L2",
                     "MCC sum [Informative]"="M1",
                     "MCC median [Informative]"="M2",
                     "Confidence sum [Not recommended]"="C1",
                     "Confidence median [Not recommended]"="C2"), width = '500px'),
       br(),
       h3("Here is the list of suggested substructures:"), 
       br(),
       dataTableOutput("table2"),
       br(),
       downloadButton("annotated_substructures", "Download substructures",style='padding:6px; font-size:150%')),
       
       column(5,
       br(),
       br(),
       br(),
       br(),
       plotOutput("plot_selected", width = '750px', height = "900px"),
       br(), offset=1)),
    
   tabPanel("Help",includeMarkdown("Help.Rmd")),
   tabPanel("About",includeMarkdown("About.Rmd"))
))