# Run packages



library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinybusy)
library(vroom)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(DOSE)
library(enrichplot)
library(DT)
library(utils)
library(base)
library(stats)
library(S4Vectors)
library(pathview)
library(PPInfer)
library(ggupset)
library(GOSemSim)
library(AnnotationDbi)



# Adding logos with links along with the title in the header
title <- tags$p("FunIntegreR",tags$a(href='https://www.uoc.edu/portal/es/index.html',
                                     
                                     tags$img(src="uoc.jpg", height = '50', width = '50')),
                tags$a(href='http://www.idisba.es/cat',
                       tags$img(src="idisba.png", height = '50', width = '50')))

ui <- dashboardPage(skin = "blue",
                    dashboardHeader(title = title, titleWidth = 400), 
                    dashboardSidebar(disable = TRUE),
                    dashboardBody(
                        
                    
                        id = "all",
                        fluidRow(
                          
                            #create info box
                            box (title= "Information", tags$p("FunIntegreR is an application developed by Genomics and Bioinformatics Platform of IDISBA, jointly with the UOC-Universitat Oberta de Catalunya and integrates
                                                             the main types of genomic functional enrichment analysis allowing the researcher to obtain genomic results in a simpler and more interactive way"),
                                                          
                                 status ="success",solidHeader = TRUE, collapsible = TRUE, width=12
                                 ),
                            #create input data box
                            
                            box(
                                title = "Input Data", width = 3, background = "light-blue", height=100, "Follow the next steps to perform the enrichment functional gene analysis" 
                            #create input results box    
                            ),
                            
                            box(
                                title = "Results", width = 9, background = "orange", height=100, "Visualize the results of your analysis" 
                                    
                                
                                
                            )
                              
                           
                                  ),
                        
                        
                        fluidRow(
                         
                        #Start Input Data function   
                            ##Gene file input box
                          
                            column(width = 3,
                                   box(
                                       title = " 1. Introduce Gene List", width = NULL, status = "primary",solidHeader = TRUE , 
                                       
                                      
                                       fileInput(inputId="file",label= "Introduce a file (.csv) with ID Genes and fold-change values (without headers)",buttonLabel = "Browse...",placeholder = "No file selected",
                                                 accept = c(".csv", ".txt"))
                                       
                                       
                                     
                                     ),
                                   
    
                                   ##Select analysis box
                                   
                                   box(
                                     title = " 2. Select Analysis", width = NULL, status = "primary",solidHeader = TRUE , collapsible= TRUE, collapsed=TRUE, 
                                     
                                     
                                     tabBox( width=12,
                                             
                                             
                                             id = "analysis",
                                      
                                             tabPanel(title="GO",
                                                      
                                                 radioButtons("GO_analysis", label=NULL,
                                                   choices = list("GO Over-Representation Analysis (GO-ORA)"= "GEA_GO", 
                                                                 "GO Gene Set Enrichment Analysis (GO-GSEA)"= "GSEA_GO"
                                                                 
                                                                  ),
                                                                  selected="GEA_GO")
                                             ),
                                             
                                             tabPanel(title="KEGG",
                                                      radioButtons("KEGG_analysis",label=NULL,
                                                                   choices = list(
                                                                                  "KEGG Over-Representation Analysis (KEGG-ORA)" = "KEGG_ORA", 
                                                                                  "KEGG Gene Set Enrichment Analysis (KEGG-GSEA)"= "KEGG_GSEA"
                                                                                  
                                                                   ),
                                                                   selected="KEGG_ORA") 
                                             ),
                                             
                                             tabPanel(title="Disease",
                                                      radioButtons("Disease_analysis",label=NULL,
                                                                   choices = list(
                                                                                  "Disease Ontology Enrichment Analysis (DO-EA)"="DO_EA",
                                                                                  "Disease Ontology Gene Set Enrichment AnalysisN(DO-GSEA)"="DO_GSEA"
                                                                   ),
                                                                   selected="DO_EA") 
                                             )
                                     
                                   )),
                                   
                                   ##Set parameters box
                                   
                                   box(
                                     title = " 3. Set the analysis parameters", width = NULL, status = "primary",solidHeader = TRUE , collapsible= TRUE, collapsed=TRUE,
                                     
                                     
                                     tabBox( width=12,
                                       
                                       
                                       id = "parameters",
                                       tabPanel(title="GO", 
                                                
                                                selectInput(inputId = "OrgDb",
                                                            label = "Organism",
                                                            choices = c("Homo sapiens"="org.Hs.eg.db","Mus musculus" ="org.Mm.eg.db")),
                                                
                                                
                                                selectInput(inputId = "ont",
                                                            label = "Subontologies(GO analysis)",
                                                            choices = c("MF", "BP", "CC"),
                                                            selected="BP"),
                                                
                                                
                                                selectInput(inputId = "keyType",
                                                            label = "Keytype of input gene",
                                                            choices = c(keytypes(org.Hs.eg.db), "MGI", "COMMON", "DESCRIPTION",
                                                                        "INTERPRO", "ORF", "SGD", "SMART"),
                                                            selected = "SYMBOL"),
                                                numericInput(inputId = "nPerm",
                                                             label = "Permutation numbers",
                                                             value = 1000,
                                                             step = 100),
                                                numericInput(inputId = "minGSSize_gse",
                                                             label = "Minimal size of each geneSet for analyzing",
                                                             value = 100,
                                                             step = 10),
                                                numericInput(inputId = "maxGSSize_gse",
                                                             label = "Maximal size of genes annotated for testing",
                                                             value = 500,
                                                             step = 1),
                                                numericInput(inputId = "pvalueCutoff",
                                                             label = "p-value cutoff",
                                                             value = 0.05,
                                                             step = 0.01)
                                                
                                               
                                                
                                                
                                                ),
                                                
                                                 
                                       
                                       tabPanel(title="KEGG",
                                                selectInput(inputId = "orgkegg",
                                                            label = "Organism",
                                                            choices = c("Homo sapiens"="hsa", "Mus musculus"="mmu")),
                                                
                                                selectInput(inputId = "keykegg",
                                                            label = "KEGG keyType",
                                                            choices = c("kegg", "uniprot","ncbi-geneid"), selected="uniprot"),
                                                
                                                numericInput(inputId = "nPerm_kegg",
                                                             label = "Permutation numbers",
                                                             value = 1000,
                                                             step = 100),
                                                
                                                numericInput(inputId = "minGSSize_gse_kegg",
                                                             label = "Minimal size of each geneSet for analyzing",
                                                             value = 10,
                                                             step = 10),
                                                
                                                numericInput(inputId = "maxGSSize_gse_kegg",
                                                             label = "Maximal size of genes annotated for testing",
                                                             value = 500,
                                                             step = 1),
                                                numericInput(inputId = "pvalueCutoff_kegg",
                                                             label = "p-value cutoff",
                                                             value = 0.05,
                                                             step = 0.01)
                                              
                                                
                                                  ),
                                       
                                       tabPanel(title="Disease", 
                                                
                                                selectInput(inputId = "keyType_do",
                                                            label = "Keytype of input gene",
                                                            choices = c(keytypes(org.Hs.eg.db), "MGI", "COMMON", "DESCRIPTION",
                                                                        "INTERPRO", "ORF", "SGD", "SMART"),
                                                            selected = "SYMBOL"),
                                 
                                                
                                                
                                                numericInput(inputId = "nPerm_do",
                                                             label = "Permutation numbers",
                                                             value = 1000,
                                                             step = 100),
                                                numericInput(inputId = "minGSSize_do",
                                                             label = "Minimal size of each geneSet for analyzing",
                                                             value = 5,
                                                             step = 10),
                                                numericInput(inputId = "maxGSSize_do",
                                                             label = "Maximal size of genes annotated for testing",
                                                             value = 500,
                                                             step = 1),
                                                numericInput(inputId = "pvalueCutoff_do",
                                                             label = "p-value cutoff",
                                                             value = 0.05,
                                                             step = 0.01)
                                                
                                                
                                                
                                                )
                                       
                                          
                                       
                                       
                                       
                                    
                                        )
                                     
                                     
                                    
                                    
                                   ),
                                   
                                   ##Submit box
                                     
                                   box(
                                     title = " 4. Submit Analysis", width = NULL, status = "primary",solidHeader = TRUE , collapsible= TRUE, collapsed=TRUE,
                                     
                                     tabBox( width=12,
                                             
                                             
                                             id = "submit_analysis",
                                             tabPanel(title="GO",
                                                       
                                            
                                                      
                                                      actionButton("goButton_GO_GEA", "Submit GO-GEA", class = "btn-primary"),
                                                      actionButton("goButton_GO_GSEA", "Submit GO-GSEA", class = "btn-primary")
                                                      
                                             ),
                                             
                                             tabPanel(title="KEGG",
                                                      actionButton("goButton_KEGG_GEA", "Submit KEGG-GEA", class = "btn-primary"),
                                                      actionButton("goButton_KEGG_GSEA", "Submit KEGG-GSEA", class = "btn-primary")
                                             ),
                                             
                                             
                                             tabPanel(title="Disease",
                                                      actionButton("goButton_DO_GEA", "Submit DO-EA", class = "btn-primary"),
                                                      actionButton("goButton_DO_GSEA", "Submit DO-GSEA", class = "btn-primary")
                                             )
                                             
                                             
                                            
                                   )
                                   
                                   ),
                                   
                                   ##Visualitzations box
                                   
                                   box(
                                     title = " 5. Select Visualizations", width = NULL, status = "primary",solidHeader = TRUE , collapsible= TRUE, collapsed=TRUE,
                                  
                                     checkboxInput("table", "Analysis Table", TRUE ),
                                     checkboxInput("dotplot_v", "Dot plot", FALSE ),
                                     checkboxInput("barplot_v","Bar Plot", FALSE),
                                     checkboxInput("enrichplot_v","Enrichment map", FALSE),
                                     checkboxInput("heatplot_v","Heatplot", FALSE),
                                     checkboxInput("network_v","Network", FALSE),
                                     checkboxInput("upsetplot_v","Upset Plot", FALSE)
                                                        
                                   ),
                                     
                                   ##Pathway box
                                     
                                   box(
                                     title = "6. KEGG Pathway Visualization", width = NULL, status = "primary",solidHeader = TRUE , collapsible= TRUE, collapsed=TRUE,
                                     
                                     textInput(inputId = "pathwayId",
                                               label = "The KEGG pathway ID"),
                                     textInput(inputId = "species",
                                               label = "Target species",
                                               value = "hsa"),
                                     
                                     shiny::actionButton(inputId = "kegg_id_path", label = "Open KEGG Pathway",  class = "btn-primary" )
                                     
                                     
                                   )
                                   
                                   
                                
                                   )  
                                     
                                  
                                   
                                   ,
                            
                            #Start Results part
                            column(width = 9,
                                   
                                  
                                   
                                   #create analysis table box
                                   box(
                                       title = "Analysis Table", width = NULL, status = "warning",solidHeader = TRUE, collapsible = TRUE,
                                       withSpinner(DT::dataTableOutput(outputId = "resultable")) 
                                       
                                   ),
                                   
                                   
                                   #create dotplot box
                                   
                                   
                                   box(
                                       
                                       title = "Dot Plot",
                                       width = NULL, solidHeader = TRUE,status = "warning", collapsible = TRUE,collapsed=TRUE,
                                       withSpinner(plotOutput("dotplot")),
                                       
                                       downloadButton(outputId = "download_dotplot", label="Download")
                                       
                                      ),
                                   
                                   #create barplot box
                                   box(
                                     title = "Bar plot", width = NULL, status = "warning",solidHeader = TRUE,  collapsible = TRUE,collapsed=TRUE,
                                     withSpinner(plotOutput(outputId = "barplot")),
                                     downloadButton(outputId = "download_barplot", label="Download")
                                   ),
                                   #create enrichment map box
                                   
                                   
                                   box(
                                     
                                     title = "Enrichment map",
                                     width = NULL, solidHeader = TRUE,status = "warning", collapsible = TRUE,collapsed=TRUE,
                                     withSpinner(plotOutput("enrichmap")),
                                     
                                     downloadButton(outputId = "download_enrichmap", label="Download")
                                     
                                   ),
                                   #create heatplot box
                                   box(
                                     
                                     title = "Heatplot",
                                     width = NULL, solidHeader = TRUE,status = "warning", collapsible = TRUE,collapsed=TRUE,
                                     withSpinner(plotOutput("heatplot")),
                                     
                                     downloadButton(outputId = "download_heatplot", label="Download")
                                     
                                   ),
                                   #create network box
                                   box(
                                     
                                     title = "Network",
                                     width = NULL, solidHeader = TRUE,status = "warning", collapsible = TRUE,collapsed=TRUE,
                                     withSpinner(plotOutput("network")),
                                     
                                     downloadButton(outputId = "download_network", label="Download")
                                     
                                   ),
                                   #create upsetplot box
                                   box(
                                     
                                     title = "Upset Plot",
                                     width = NULL, solidHeader = TRUE,status = "warning", collapsible = TRUE,collapsed=TRUE,
                                     withSpinner(plotOutput("upsetplot")),
                                     
                                     downloadButton(outputId = "download_upsetplot", label="Download")
                                     
                                   )

                                     
                                   
                                     
                                                                   
                            
                                )
                                   
                            
                      )
                        
                )
                    
                    
                    
)





                         


                            

                                   
                                   
                 

