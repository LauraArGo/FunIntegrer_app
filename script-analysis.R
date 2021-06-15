
#install packages


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("GO.db")
BiocManager::install("DO.db")
BiocManager::install("DOSE")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("PPInfer")
BiocManager::install("GOSemSim")
install.packages("ggplot2")
install.packages("ggupset")
install.packages("phylobase")



library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(DOSE)
library(enrichplot)
library(pathview)
library(PPInfer)
library(ggupset)
library(phylobase)
library(GOSemSim)
library(AnnotationDbi)

#The geneList contains three features:
  
#1. numeric vector: fold change or other type of numerical variable
#2. named vector: every number was named by the corresponding gene ID
#3.sorted vector: number should be sorted in decreasing order

#Suppose you are importing your own data from a csv file and the file contains two columns, 
#one for gene ID (no duplicated allowed) and another one for fold change, you can prepare your
#own geneList via the following command:


data <- read.table("datos-prueba.csv",header = FALSE, sep=c(";"))
## assume that 1st column is ID
## 2nd column is fold change

## feature 1: numeric vector
geneList = data[,2]

## feature 2: named vector
names(geneList) = data[,1]


## feature 3: sort vector

geneList = sort(geneList, decreasing = TRUE)

head(geneList)



#Go analysis


genes<-as.character(data[,1])

ego<-enrichGO(genes, OrgDb='org.Hs.eg.db', keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)




#barplot

barplot(ego, showCategory=20)

#dotplot
dotplot(ego, showCategory=20)

#emapplot
d <- godata('org.Hs.eg.db', ont="BP")
x<-ego
ego2 <- pairwise_termsim(x)
emapplot(ego2, showCategory=20)
emapplot_cluster(ego2,showCategory=20)


#heatplot



heatplot(ego, foldChange=geneList)

#Network
cnetplot(ego, foldChange=geneList)

#upsetplot
upsetplot(ego)


##Gene Set Enrichment Analysis (GSEA) of GO



ego3 <- gseGO(geneList= geneList,OrgDb = get("org.Hs.eg.db"), keyType = "SYMBOL", ont = "BP", nPerm= 1000,minGSSize= 100,
              maxGSSize = 500,  pvalueCutoff = 0.05, verbose= FALSE)

upsetplot(ego3)
dotplot(ego3)
ego4<-pairwise_termsim(ego3)
emapplot(ego4)
cnetplot(ego3, foldChange=geneList)

e_GO<-function(geneList){
  
  
  
  gseGO(geneList= x,OrgDb = get("org.Hs.eg.db"), keyType = "SYMBOL", ont = "BP", nPerm= 1000,minGSSize= 100,
        maxGSSize = 500,  pvalueCutoff = 0.05, verbose= FALSE)

 
    return(gseGO)
}  


e_go2<-as.data.frame(e_GO()@result)

ego3_table<-as.data.frame(ego3)



GSEA.barplot(ego3_table, category = "Description", score = "p.adjust",pvalue = "setSize",top = 20, sort  = "p.adjust")



## KEGG analysis

#KEGG ORA

#convert gene symbol to UNIPROT

gene.df <- bitr(genes, fromType = "SYMBOL",
                toType = "UNIPROT",
                OrgDb = org.Hs.eg.db)
head(gene.df)

#KEGG enrichment
e_KEGG<-enrichKEGG(gene.df[,2], organism = "hsa", keyType = "uniprot", pvalueCutoff = 0.05)

head(e_KEGG)

dotplot(e_KEGG)
barplot(e_KEGG)
upsetplot(e_KEGG)

e_KEGG2 <- pairwise_termsim(e_KEGG)
emapplot(e_KEGG2, showCategory=20)

##KEGG Gene Set Enrichment Analysis

#creat GeneList (keytype:UNIPROT) with Fold-change

gene.df2<-merge(gene.df, data, by.x="SYMBOL", by.y=c("V1"))

geneList_kk = gene.df2[,3]

names(geneList_kk) = gene.df2[,2]

geneList_kk = sort(geneList_kk, decreasing = TRUE)

head(geneList_kk)

#KEGG GSEA function

gsea_KEGG <- gseMKEGG(geneList = geneList_kk, 
                 organism = 'hsa', keyType = "uniprot",nPerm = 1000, minGSSize = 10, maxGSSize = 500,
                 pvalueCutoff = 0.05)

head(gsea_KEGG)
dotplot(gsea_KEGG)

## Visualize enriched KEGG pathways
browseKEGG(e_KEGG, 'hsa01240')


hsa04010 <- pathview(gene.data  = geneList_kk,
                     pathway.id = "hsa01240",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList_kk)), cpd=1))

#Disease analysis


##convert gene symbol to ENTREZ ID

gene.do <- bitr(genes, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

##creat GeneList (keytype:ENTREZ ID) with Fold-change

gene.do2<-merge(gene.do, data, by.x="SYMBOL", by.y=c("V1"))

geneList_do = gene.do2[,3]

names(geneList_do) = gene.do2[,2]

geneList_do = sort(geneList_do, decreasing = TRUE)

head(geneList_do)

##enrichdo function,  Disease Ontology Enrichment Analysis (DO-EA) 

e_DO <- enrichDO(gene          = gene.do[,2],
              ont           = "DO",
              pvalueCutoff  = 0.05,
              minGSSize     = 5,
              maxGSSize     = 500)
head(e_DO)

dotplot(e_DO)

#Disease Ontology Gene Set Enrichment Analysis (DO-GSEA)


gsea_DO <- gseDO(geneList_do,
           minGSSize     = 120,
           pvalueCutoff  = 0.2,
           pAdjustMethod = "BH",
           verbose       = FALSE)
head(gsea_DO, 3)

dotplot(gsea_DO)

