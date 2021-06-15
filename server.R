

library(shiny)
library(shinydashboard)


server <- function(input, output, session){

#1. INPUT DATA 
  
## Load data from genes file. assume that 1st column is ID. & 2nd column is fold change
  
    gene_data <- reactive({
      validate(
        need(input$file !="", "Please load the gene file")
      )
      
      gene_data1 = input$file
      data1 = read.table(gene_data1$datapath, sep=";")
      return(data1)
    
    
  })
    
    
    
## Create genesList with 2 vectors: IDs genes vector  and numeric vector (fold-change)
   
    genesList <- reactive({

      
      ## feature 1: numeric vector
      
      genesList = gene_data()$V2
      
      ## feature 2: named vector
      names(genesList) = gene_data()$V1
      
      
      ## feature 3: sort vector
      
      genesList = sort(genesList, decreasing = TRUE)
      
      return(genesList)
      
    })
    
 

    
#2. ANALYSIS
    
  ##  EnrichGO function
  e_GO <- eventReactive(input$goButton_GO_GEA,{
    
    
    table<-as.character(gene_data()$V1)
    
     e_GO <- enrichGO(
      gene = table, OrgDb = get(input$OrgDb), keyType = input$keyType, ont = input$ont,
        pvalueCutoff = as.numeric(input$pvalueCutoff)
      )
    
  })
  
  ##  GO GSEA function
  gsea_GO <- eventReactive(input$goButton_GO_GSEA,{
    
   
    
    table<-genesList()
    gsea_GO <- gseGO(
      gene = table, OrgDb = get(input$OrgDb), keyType = input$keyType, ont = input$ont, nPerm= input$nPerm, minGSSize= input$minGSSize_gse,
              maxGSSize = input$maxGSSize_gse, pvalueCutoff = as.numeric(input$pvalueCutoff)
      )
    
  })
  

  
  #KEGG enrichment function
  
  
  e_KEGG <- eventReactive(input$goButton_KEGG_GEA,{
    
    
    genes<-gene_data()$V1
    data <- bitr(genes, fromType = "SYMBOL",
                    toType = "UNIPROT",
                    OrgDb = org.Hs.eg.db)
    data<-data[,2]
    
    e_KEGG<-enrichKEGG(data, organism = input$orgkegg, keyType = input$keykegg, pvalueCutoff = input$pvalueCutoff_kegg)
  
    
    
    })
  
  ## KEGG GSEA function
  
geneList_kk<- reactive({
    
    table<-gene_data()
    
    genes<-as.character(table[,1])
    
    genes <- bitr(genes, fromType = "SYMBOL",
                  toType = "UNIPROT",
                  OrgDb = org.Hs.eg.db)
    
    genes_kk<-merge(genes, table, by.x="SYMBOL", by.y=c("V1"))
    
    geneList_kk = genes_kk[,3]
    
    names(geneList_kk) = genes_kk[,2]
    
    geneList_kk = sort(geneList_kk, decreasing = TRUE)
    
    return(geneList_kk)
    
})



  
  gsea_KEGG <- eventReactive(input$goButton_KEGG_GSEA,{
    
    
     
      table <- geneList_kk()
      gsea_KEGG <- gseMKEGG(geneList = table,
                          organism = input$orgkegg , keyType = input$keykegg, nPerm = input$nPerm_kegg, minGSSize = input$minGSSize_gse_kegg, maxGSSize = input$maxGSSize_gse_kegg,
                          pvalueCutoff = input$pvalueCutoff_kegg)
    
  })
  
  ##  Disease Analysis
  
    ###DO enrichment function
  
  
   e_DO<- eventReactive(input$goButton_DO_GEA,{
    
  
      genes<-gene_data()$V1
      data <- bitr(genes, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
      data<-data[,2]
      
      e_DO<-enrichDO(data, ont           = "DO",
                    
                     minGSSize     = input$minGSSize_do,
                     maxGSSize     = input$maxGSSize_do,
                     pvalueCutoff = input$pvalueCutoff_do)
    
    
    
  })
  
  #DO-GSEA function
    
  geneList_do<- reactive({
    
    
    
    if(input$keyType != "ENTREZID"){
      
    table<-gene_data()
      
    genes<-as.character(table[,1])
    
    genes <- bitr(genes, fromType = input$keyType_do,
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
  
    
    genes_do<-merge(genes, table, by.x=input$keyType_do, by.y=c("V1"))
    geneList_do = genes_do[,3]
    
    names(geneList_do) = genes_do[,2]
    
    geneList_do = sort(geneList_do, decreasing = TRUE)
    
    
    }
    
    return(geneList_do)
    
  })
  
  gsea_do <- eventReactive(input$goButton_DO_GSEA,{
    
  
      
      table <- geneList_do()
      gsea_do <- gseDO(table,
                       minGSSize     = 120,
                       pvalueCutoff  = 0.2,
                       pAdjustMethod = "BH",
                       verbose       = FALSE)
     
  })
  
  
  #3. PARAMETERS 
  
   ## Set Parameters Tab Panel depending on analysis selected:
  observe({
    
    a<- input$analysis
    
    if (a == "GO" ) {
      updateTabsetPanel(session, "parameters", selected = "GO")
    } else if (a == "KEGG" ){
      updateTabsetPanel(session, "parameters", selected = "KEGG")
    } else if (a == "Disease"){
      updateTabsetPanel(session, "parameters", selected = "Disease")
    }
    
    
    
  })

# 4. SUBMIT ANALYSIS  
  
  ## Set SUBMIT Tab Panel depending on parameters selected:
  
  observe({
    
    a<- input$parameters
    
    if (a == "GO" ) {
      updateTabsetPanel(session, "submit_analysis", selected = "GO")
    } else if (a == "KEGG" ){
      updateTabsetPanel(session, "submit_analysis", selected = "KEGG")
    } else if (a == "Disease"){
      updateTabsetPanel(session, "submit_analysis", selected = "Disease")
    }
    
  })
    
  
# 5. VISUALIZATIONS
  
  
  ##Output table
  
 
  output$resultable <- DT::renderDataTable({
     
     req(input$table)
     
     a<- input$analysis
     
    if(a == "GO" && input$GO_analysis== "GEA_GO"){
        my_data<-e_GO()@result
        DT::datatable(my_data,extensions = 'Buttons',
                  options = list(
                    paging = TRUE,
                    searching = TRUE,
                    fixedColumns = TRUE,
                    autoWidth = TRUE,
                    scrollX = TRUE,
                    ordering = TRUE,
                    dom = 'Bfrtip',
                    buttons = c('copy', 'csv', 'excel')
                  ))
       }
     else if(a == "GO" && input$GO_analysis== "GSEA_GO"){
        my_data<-gsea_GO()@result
        DT::datatable(my_data,extensions = 'Buttons',
                  options = list(
                    paging = TRUE,
                    searching = TRUE,
                    fixedColumns = TRUE,
                    autoWidth = TRUE,
                    scrollX = TRUE,
                    ordering = TRUE,
                    dom = 'Bfrtip',
                    buttons = c('copy', 'csv', 'excel')
                  ))
      }
     else if(a == "KEGG" && input$KEGG_analysis== "KEGG_ORA"){
      my_data<-e_KEGG()@result
      DT::datatable(my_data,extensions = 'Buttons',
                    options = list(
                      paging = TRUE,
                      searching = TRUE,
                      fixedColumns = TRUE,
                      autoWidth = TRUE,
                      scrollX = TRUE,
                      ordering = TRUE,
                      dom = 'Bfrtip',
                      buttons = c('copy', 'csv', 'excel')
                    ))
      }
    
     else if(a == "KEGG" && input$KEGG_analysis== "KEGG_GSEA"){
      my_data<- gsea_KEGG()@result
      DT::datatable(my_data,extensions = 'Buttons',
                    options = list(
                      paging = TRUE,
                      searching = TRUE,
                      fixedColumns = TRUE,
                      autoWidth = TRUE,
                      scrollX = TRUE,
                      ordering = TRUE,
                      dom = 'Bfrtip',
                      buttons = c('copy', 'csv', 'excel')
                    ))
     }
     
      else if(a == "Disease" && input$Disease_analysis== "DO_EA"){
        my_data<- e_DO()@result
        DT::datatable(my_data,extensions = 'Buttons',
                      options = list(
                        paging = TRUE,
                        searching = TRUE,
                        fixedColumns = TRUE,
                        autoWidth = TRUE,
                        scrollX = TRUE,
                        ordering = TRUE,
                        dom = 'Bfrtip',
                        buttons = c('copy', 'csv', 'excel')
                      ))
    }
     else if(a == "Disease" && input$Disease_analysis== "DO_GSEA"){
       my_data<- gsea_do()@result
       DT::datatable(my_data,extensions = 'Buttons',
                     options = list(
                       paging = TRUE,
                       searching = TRUE,
                       fixedColumns = TRUE,
                       autoWidth = TRUE,
                       scrollX = TRUE,
                       ordering = TRUE,
                       dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel')
                     ))
     }
    
    })
  
 
  

  
  ##Output dotplot
   
output$dotplot <- renderPlot({ 
    
    req(input$dotplot_v)
    
    a<- input$analysis
    
    if ( a== "GO" && "GEA_GO" %in% input$GO_analysis) {
      
        dotplot(e_GO(), showCategory=20)
      
    }
    
    else if (a== "GO" && "GSEA_GO" %in% input$GO_analysis) {
      
              dotplot(gsea_GO(), showCategory=20)
      
      
    }
    
    else if (a== "KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
      
             dotplot(e_KEGG(), showCategory=20)
      
      
    }
    
    else if (a== "KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
      
              dotplot(gsea_KEGG(), showCategory=20)
      
      
    }
    
    else if (a== "Disease" && "DO_EA" %in% input$Disease_analysis) {
      
      dotplot(e_DO(), showCategory=20)
      
      
    }
    
    else if (a== "Disease" && "DO_GSEA" %in% input$Disease_analysis) {
      
      dotplot(gsea_do(), showCategory=20)
      
      
    }
    
  })
  
##Output barplot


### For GSEA.barplot, it's necesary create a function return GSEA results in a data.frame format


#### GO.GSEA dataframe

gsea_GO2<-reactive({
  
  table<-genesList()
  gsea_GO2 <- gseGO(
    gene = table, OrgDb = get(input$OrgDb), keyType = input$keyType, ont = input$ont, nPerm= input$nPerm, minGSSize= input$minGSSize_gse,
    maxGSSize = input$maxGSSize_gse, pvalueCutoff = as.numeric(input$pvalueCutoff)
  )
  gsea_GO2<-as.data.frame(gsea_GO2)
  
})
#### KEGG.GSEA dataframe
gsea_KEGG2 <- reactive({
  
  
  
  table <- geneList_kk()
  gsea_KEGG2 <- gseMKEGG(geneList = table,
                         organism = input$orgkegg , keyType = input$keykegg, nPerm = input$nPerm_kegg, minGSSize = input$minGSSize_gse_kegg, maxGSSize = input$maxGSSize_gse_kegg,
                         pvalueCutoff = input$pvalueCutoff_kegg)
  gsea_KEGG2<-as.data.frame(gsea_KEGG2)
  
})

#### DO.GSEA dataframe


gsea_do2 <- reactive({
  
  
  
  table <- geneList_do()
  gsea_do2 <- gseDO(table,
                    minGSSize     = 120,
                    pvalueCutoff  = 0.2,
                    pAdjustMethod = "BH",
                    verbose       = FALSE)
  gsea_do2<-as.data.frame(gsea_do2)
  
})


####output barplot
output$barplot <- renderPlot({
  
  req(input$barplot_v)
  
  a<- input$analysis
  table1<-gsea_GO2()
  table2<-gsea_KEGG2()
  table3<-gsea_do2()
  
  
  
  if (a=="GO" && "GEA_GO" %in% input$GO_analysis) {
    
    barplot(e_GO(), showCategory=20)
    
  }
  
  else if (a=="GO" && "GSEA_GO" %in% input$GO_analysis) {
    
    
    
    GSEA.barplot(table1, category = "Description", score = "NES",pvalue = "p.adjust",top = 20, sort  = "NES")
    
    
  }
  
  else if (a=="KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
    
    barplot(e_KEGG(), showCategory=20)
    
    
  }
  
  else if (a=="KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
    
    
    
    GSEA.barplot(table2, category = "Description", score = "NES",pvalue = "p.adjust",top = 20, sort  = "NES")
    
    
  }
  
  else if (a=="Disease" && "DO_EA" %in% input$Disease_analysis) {
    
    barplot(e_DO(), showCategory=20)
    
    
  }
  
  else if (a=="Disease" && "DO_GSEA" %in% input$Disease_analysis) {
    
    
    
    GSEA.barplot(table3, category = "Description", score = "NES",pvalue = "p.adjust",top = 20, sort  = "NES")
    
    
  }
  
  
})

e_GO3 <- reactive({
  
  
  table<-as.character(gene_data()$V1)
  
  e_GO3 <- enrichGO(
    gene = table, OrgDb = get(input$OrgDb), keyType = input$keyType, ont = input$ont,
    pvalueCutoff = as.numeric(input$pvalueCutoff)
  )
  
})
  
  ##Output enrichmap


       
  output$enrichmap <- renderPlot({ 
      
      req(input$enrichplot_v)
      
      a<- input$analysis
     
      
      d1<-e_GO()
      data1 <- enrichplot::pairwise_termsim(d1)
      gsea_GO3<-pairwise_termsim(gsea_GO())
      e_KEGG3 <- pairwise_termsim(e_KEGG())
      gsea_KEGG3<- pairwise_termsim(gsea_KEGG())
      e_DO3<- pairwise_termsim(e_DO())
      gsea_do3<-pairwise_termsim(gsea_do())
      
      
      if (a=="GO" && "GEA_GO" %in% input$GO_analysis) {
        
        
        
        enrichplot::emapplot(data1, showCategory=20)
        
      }
      
      else if (a=="GO" && "GSEA_GO" %in% input$GO_analysis) {
        
        emapplot(gsea_GO2, showCategory=20)
        
        
      }
      
      else if (a=="KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
        
        emapplot(e_KEGG3, showCategory=20)
        
        
      }
      
      else if (a=="KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
        
        emapplot(gsea_KEGG3, showCategory=20)
        
        
      }
      else if (a=="Disease" && "DO_EA" %in% input$Disease_analysis) {
        
        emapplot(e_DO3, showCategory=20)
        
        
      }
      
      else if (a=="Disease" && "DO_GSEA" %in% input$Disease_analysis) {
        
        emapplot(gsea_do3, showCategory=20)
        
        
      }
  })
    

        
    
   
    ##Output heatplot
    
      
    
  output$heatplot<- renderPlot({ 
      
      req(input$heatplot_v)
      
      a<- input$analysis
      
      if (a== "GO" && "GEA_GO" %in% input$GO_analysis) {
        
        table_fold<-genesList()
        
        heatplot(e_GO(), foldChange=table_fold, showCategory=20)
        
      }
      
      else if (a== "GO" && "GSEA_GO" %in% input$GO_analysis) {
        
        table_fold<-genesList()
        
        heatplot(gsea_GO(), foldChange=table_fold, showCategory=20)
        
        
      }
      
      else if (a== "KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
        
        table_fold<-geneList_kk()
        
        heatplot(e_KEGG(), foldChange=table_fold, showCategory=20)
        
        
      }
      
      else if (a== "KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
        
        table_fold<-geneList_kk()
        
        heatplot(gsea_KEGG(), foldChange=table_fold, showCategory=20)
        
        
        
      }
      
      else if (a== "Disease" && "DO_EA" %in% input$Disease_analysis) {
        
        table_fold<-geneList_do()
        
        heatplot(e_DO(), foldChange=table_fold, showCategory=20)
        
        
      }
      
      else if (a== "Disease" && "DO_GSEA"%in% input$Disease_analysis) {
        
        table_fold<-geneList_do()
        
        heatplot(gsea_do(), foldChange=table_fold, showCategory=20)
        
        
        
      }
      
    })
    
    
    ##Output network
  output$network<- renderPlot({ 
      
       req(input$network_v)
      
       a<- input$analysis
      
        
        if (a== "GO" && "GEA_GO" %in% input$GO_analysis) {
          table_fold<-genesList()
          
          cnetplot(e_GO(), foldChange=table_fold, showCategory=20)
          
        }
        
        else if (a== "GO" && "GSEA_GO" %in% input$GO_analysis) {
          
          table_fold<-genesList()
          
          cnetplot(gsea_GO(), foldChange=table_fold, showCategory=20)
          
          
        }
        
        else if (a== "KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
          
          table_fold<-geneList_kk()
          
          cnetplot(e_KEGG(), foldChange=table_fold, showCategory=20)
          
          
        }
        
        else if (a== "KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
          
          table_fold<-geneList_kk()
          
          cnetplot(gsea_KEGG(), foldChange=table_fold, showCategory=20)
          
          
          
        }
      
      else if (a== "Disease" && "DO_EA" %in% input$Disease_analysis) {
        
        table_fold<-geneList_do()
        
        cnetplot(e_DO(), foldChange=table_fold, showCategory=20)
        
        
      }
      
      else if (a== "Disease" && "DO_GSEA" %in% input$Disease_analysis) {
        
        table_fold<-geneList_do()
        
        cnetplot(gsea_do(), foldChange=table_fold, showCategory=20)
        
        
        
      }
        
      
    })
    

   ##Output upset plot
  output$upsetplot<- renderPlot({ 
     
     
     req(input$upsetplot_v)
     
     a<- input$analysis

     
     if (a=="GO" && "GEA_GO" %in% input$GO_analysis) {
       
       upsetplot(e_GO())
       
     }
     
     else if (a=="GO" && "GSEA_GO" %in% input$GO_analysis) {
       
       upsetplot(gsea_GO())
       
       
     }
     
     else if (a=="KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
       
       upsetplot(e_KEGG())
       
       
     }
     
     else if (a=="KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
       
       upsetplot(gsea_KEGG())
       
       
     }
     
     else if (a=="Disease" && "DO_EA" %in% input$Disease_analysis) {
       
       upsetplot(e_DO())
       
       
     }
     
     else if (a=="Disease" && "DO_GSEA" %in% input$Disease_analysis) {
       
       upsetplot(gsea_do())
       
     }
  
    })

   ##Output kegg pathway
 
  
  observeEvent(input$kegg_id_path,{
    id<-input$pathwayId;
    
    browseURL(paste0("https://www.genome.jp/pathway/",id))
  }) 

   
   
# DOWNLOADS
    
    ## Download dotplot  
    output$download_dotplot <- downloadHandler(
      filename = ("dotplot.png"),
      content = function(file) {
        
        req(input$dotplot_v)
        a<- input$analysis 
      
        if (a=="GO" && "GEA_GO" %in% input$GO_analysis) {
          ggsave(file, clusterProfiler::dotplot(e_GO(),showCategory=20),width=14, height=8)
        }
        
         else if (a=="GO" && "GSEA_GO" %in% input$GO_analysis) {
          ggsave(file, clusterProfiler::dotplot(gsea_GO(), showCategory=20),width=14, height=8)
         }
        else if (a=="KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
          ggsave(file, clusterProfiler::dotplot(e_KEGG(), showCategory=20),width=14, height=8)
        }
        else if (a=="KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
          ggsave(file, clusterProfiler::dotplot(gsea_KEGG(), showCategory=20),width=14, height=8)
        }
        
        else if (a=="Disease" && "DO_EA" %in% input$Disease_analysis) {
          ggsave(file, clusterProfiler::dotplot(e_DO(), showCategory=20),width=14, height=8)
        }
        else if (a=="Disease" && "DO_GSEA" %in% input$Disease_analysis) {
          ggsave(file, clusterProfiler::dotplot(gsea_do(), showCategory=20),width=14, height=8)
        }
      }
    )
   
 
       ## Download barplot  

    
    output$download_barplot <- downloadHandler(
      filename = ("barplot.png"),
      content = function(file) {
        
        req(input$barplot_v)
        a<- input$analysis 
        
        table1<-as.data.frame(gsea_go())
        table2<-as.data.frame(gsea_KEGG())
        table3<-as.data.frame(gsea_do())
        
        
        
        if (a=="GO" && "GEA_GO" %in% input$GO_analysis) {
          ggsave(file, clusterProfiler::barplot(e_GO(),showCategory=20),width=14, height=8)
        }
        
        else if (a=="GO" && "GSEA_GO" %in% input$GO_analysis) {
          ggsave(file, clusterProfiler::GSEA.barplot(table1, category = "Description", score = "setSize",pvalue = "p.adjust",top = 20, sort  = "setSize"),width=14, height=8)
        }
        else if (a=="KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
          ggsave(file, clusterProfiler::barplot(e_KEGG(), showCategory=20),width=14, height=8)
        }
        else if (a=="KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
          ggsave(file, clusterProfiler::GSEA.barplot(table2, category = "Description", score = "setSize",pvalue = "p.adjust",top = 20, sort  = "setSize"),width=14, height=8)
        }
        
        else if (a=="Disease" && "DO_EA" %in% input$Disease_analysis) {
          ggsave(file, clusterProfiler::barplot(e_DO(), showCategory=20),width=14, height=8)
        }
        else if (a=="Disease" && "DO_GSEA" %in% input$Disease_analysis) {
          ggsave(file, clusterProfiler::GSEA.barplot(table3, category = "Description", score = "setSize",pvalue = "p.adjust",top = 20, sort  = "setSize"),width=14, height=8)
        }
      }
    )

    ## Download enrichplot 
    
    
    output$download_enrichmap <- downloadHandler(
      filename = ("enrichplot.png"),
      content = function(file) {
        
        req(input$enrichplot_v)
        a<- input$analysis 
        
        if (a=="GO" && "GEA_GO" %in% input$GO_analysis) {
          ggsave(file, clusterProfiler::emapplot(e_GO(),showCategory=20),width=14, height=8)
        }
        
        else if (a=="GO" && "GSEA_GO" %in% input$GO_analysis) {
          ggsave(file, clusterProfiler::emapplot(gsea_GO(), showCategory=20),width=14, height=8)
        }
        else if (a=="KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
          ggsave(file, clusterProfiler::emapplot(e_KEGG(), showCategory=20),width=14, height=8)
        }
        else if (a=="KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
          ggsave(file, clusterProfiler::emapplot(gsea_KEGG(), showCategory=20),width=14, height=8)
        }
        
        else if (a=="Disease" && "DO_EA" %in% input$Disease_analysis) {
          ggsave(file, clusterProfiler::emapplot(e_DO(), showCategory=20),width=14, height=8)
        }
        else if (a=="Disease" && "DO_GSEA" %in% input$Disease_analysis) {
          ggsave(file, clusterProfiler::emapplot(gsea_do(), showCategory=20),width=14, height=8)
        }
      }
    )
                      
    ## Download upsetplot 
    
    
    output$download_upsetplot <- downloadHandler(
      filename = ("upsetplot.png"),
      content = function(file) {
        
        req(input$upsetplot_v)
        a<- input$analysis 
        
        if (a=="GO" && "GEA_GO" %in% input$GO_analysis) {
          ggsave(file, clusterProfiler::upsetplot(e_GO(),showCategory=20),width=14, height=8)
        }
        
        else if (a=="GO" && "GSEA_GO" %in% input$GO_analysis) {
          ggsave(file, clusterProfiler::upsetplot(gsea_GO(), showCategory=20),width=14, height=8)
        }
        else if (a=="KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
          ggsave(file, clusterProfiler::upsetplot(e_KEGG(), showCategory=20),width=14, height=8)
        }
        else if (a=="KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
          ggsave(file, clusterProfiler::upsetplot(gsea_KEGG(), showCategory=20),width=14, height=8)
        }
        
        else if (a=="Disease" && "DO_EA" %in% input$Disease_analysis) {
          ggsave(file, clusterProfiler::upsetplot(e_DO(), showCategory=20),width=14, height=8)
        }
        else if (a=="Disease" && "DO_GSEA" %in% input$Disease_analysis) {
          ggsave(file, clusterProfiler::upsetplot(gsea_do(), showCategory=20),width=14, height=8)
        }
      }
    )
  
## Download heatplot    

  
  output$download_heatplot <- downloadHandler(
          filename = ("heatplot.png"),
          content = function(file) {
            
            req(input$heatplot_v)
            a<- input$analysis 
            table_fold1<-genesList()
            table_fold2<-geneList_kk()
            table_fold3<-geneList_do()
            
            if (a=="GO" && "GEA_GO" %in% input$GO_analysis) {
              ggsave(file, clusterProfiler::heatplot(e_GO(), foldChange=table_fold1,showCategory=20),width=14, height=8)
            }
            
            else if (a=="GO" && "GSEA_GO" %in% input$GO_analysis) {
              ggsave(file, clusterProfiler::heatplot(gsea_GO(), foldChange=table_fold1, showCategory=20),width=14, height=8)
            }
            else if (a=="KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
              ggsave(file, clusterProfiler::heatplot(e_KEGG(), foldChange=table_fold2, showCategory=20),width=14, height=8)
            }
            else if (a=="KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
              ggsave(file, clusterProfiler::heatplot(gsea_KEGG(), foldChange=table_fold2, showCategory=20),width=14, height=8)
            }
            
            else if (a=="Disease" && "DO_EA" %in% input$Disease_analysis) {
              ggsave(file, clusterProfiler::heatplot(e_DO(), foldChange=table_fold3, showCategory=20),width=14, height=8)
            }
            else if (a=="Disease" && "DO_GSEA" %in% input$Disease_analysis) {
              ggsave(file, clusterProfiler::heatplot(gsea_do(), foldChange=table_fold3,showCategory=20),width=14, height=8)
            }
          }
        )
        
  ## Download network    
  
  
  output$download_network <- downloadHandler(
    filename = ("network.png"),
    content = function(file) {
      
      req(input$network_v)
      a<- input$analysis 
      table_fold1<-genesList()
      table_fold2<-geneList_kk()
      table_fold3<-geneList_do()
      
      if (a=="GO" && "GEA_GO" %in% input$GO_analysis) {
        ggsave(file, clusterProfiler::cnetplot(e_GO(), foldChange=table_fold1,showCategory=20),width=14, height=8)
      }
      
      else if (a=="GO" && "GSEA_GO" %in% input$GO_analysis) {
        ggsave(file, clusterProfiler::cnetplot(gsea_GO(), foldChange=table_fold1, showCategory=20),width=14, height=8)
      }
      else if (a=="KEGG" && "KEGG_ORA" %in% input$KEGG_analysis) {
        ggsave(file, clusterProfiler::cnetplot(e_KEGG(), foldChange=table_fold2, showCategory=20),width=14, height=8)
      }
      else if (a=="KEGG" && "KEGG_GSEA" %in% input$KEGG_analysis) {
        ggsave(file, clusterProfiler::cnetplot(gsea_KEGG(), foldChange=table_fold2, showCategory=20),width=14, height=8)
      }
      
      else if (a=="Disease" && "DO_EA" %in% input$Disease_analysis) {
        ggsave(file, clusterProfiler::cnetplot(e_DO(), foldChange=table_fold3, showCategory=20),width=14, height=8)
      }
      else if (a=="Disease" && "DO_GSEA" %in% input$Disease_analysis) {
        ggsave(file, clusterProfiler::cnetplot(gsea_do(), foldChange=table_fold3,showCategory=20),width=14, height=8)
      }
    }
  )     
        
  
}       




   

    




    



