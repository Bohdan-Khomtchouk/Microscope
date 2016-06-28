# Copyright (C) 2015-2016 Bohdan Khomtchouk, James Hennessy, and Vytas Dargis-Robinson

# This file is part of MicroScope.

# ------------------------------------------------------------------------------------

library(shiny)
library(d3heatmap)
library(RColorBrewer)
library(edgeR)
library(GO.db)
library(goseq)
library(networkD3)


# frontend
ui <- shinyUI(pageWithSidebar(
  
  headerPanel("MicroScope"),
  
  sidebarPanel(
    downloadButton("downloadData", label = "Download Sample Input File"),
    fileInput("filename", "Choose File to Upload:", accept = c('.csv')),
    uiOutput("expcolumns"),
    actionButton("goButton", "Run Statistics"),
    downloadButton("downloadtable", "Download Stats Table"),
    
    conditionalPanel(
      condition = "input.goButton == true",
      selectInput("pvFDRchoose", "Choose Heatmap Statistical Criterion:", selected = "Pvalue", c("Pvalue", "FDR", "FDR and PValue" = "both")
    ),
    conditionalPanel(
      condition= "input.pvFDRchoose == 'FDR' | input.pvFDRchoose == 'both'",
      sliderInput("statFDR", "Choose Cutoff for FDR in Heatmap", min = 0.001, max = 0.1, value = 0.1)
    ),
    conditionalPanel(
      condition= "input.pvFDRchoose == 'Pvalue' | input.pvFDRchoose == 'both'",
      sliderInput("statPV", "Choose Cutoff for Pvalues in Heatmap", min = 0.001, max = 0.05, value = 0.05)
    ),
    
    checkboxInput("log2_transformed_data", "log2 Transform Heatmap"),
    selectInput("choose", "Choose Color Scheme:", c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu", "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges", "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues")),
    selectInput("dendrogram", "Apply Clustering:", c("none", "row", "column", "both")),
    numericInput("color_row_branches", "Color Row Branches:", value = 1),
    numericInput("color_column_branches", "Color Column Branches:", value = 1),
    sliderInput("xfontsize", "Choose X Font Size:", min = 0.3, max = 2, value = 0.5),
    sliderInput("yfontsize", "Choose Y Font Size:", min = 0.3, max = 2, value = 1.0),
    actionButton("goButtonHeat", "Draw Heatmap"),
    downloadButton("downloadHeatmap", "Download Heatmap")
    ),
    
    conditionalPanel(
      condition = "input.goButtonHeat == true",
      selectInput("Type", "Choose PCA Option:", selected = "Correlation Matrix", c("Covariance Matrix", "Correlation Matrix")),
      actionButton("pcaButton", "Run PCA"),
      downloadButton("downloadBiplot", "Download PCA Biplot"),
      downloadButton("downloadScreeplot", "Download PCA Screeplot")
      ),
      
    conditionalPanel(
      condition = "input.pcaButton == true",  
      selectInput("Genome", "Choose Genome Database:", selected = "mm9", c("Mouse" = "mm9", "Human" = "hg19", "Chimpanzee" = "panTro2", "Rat" = "rn4", "Worm" = "ce6", "Zebrafish" = "danRer6", "Fly" = "dm3", "Yeast" = "sacCer2", "Cow" = "bosTau4", "Dog" = "canFam2", "Anopheles gambiae" = "anoGam1", "Rhesus" = "rheMac2", "Frog" = "xenTro2", "Chicken" = "galGal3")),
      selectInput("geneRef", "Choose Gene Identifier:", selected = "geneSymbol", c("Gene Symbol" = "geneSymbol", "Ensembl ID" = "ensGene")),
      numericInput("numberGenes", "Choose How Many Top Gene Ontologies to Display:", value = 10),
      selectInput("chooseEnriched", "Stratify Top Gene Ontologies By:", c("P-value", "FDR")),
      sliderInput("cutoffP", "Choose Gene Ontology P-value Cutoff:", min = 0.001, max = 0.1, value = 0.05),
	  sliderInput("cutoffFdr", "Choose Gene Ontology FDR Cutoff:", min = 0.001, max = 0.1, value = 0.05),
	  actionButton("goData", "Do Gene Ontology Analysis")
    ),
    
    conditionalPanel(
      condition = "input.goData == true",
      actionButton("doNets", "Do Network Analysis"),
      downloadButton("downloadSimpleNetwork", "Download Network")
    )
    
  ),

  mainPanel(
    tabsetPanel(
      tabPanel("Instructions", textOutput("text1"), img(src='excel.png'), textOutput("text2"), textOutput("text3"), textOutput("text4"), textOutput("text5"), textOutput("text6"), textOutput("text7"), textOutput("text8"), textOutput("text9")),
      tabPanel("DE Analysis", tableOutput("table")),
      tabPanel("Heatmap", uiOutput("pixelation")),
      tabPanel("PCA", verbatimTextOutput("pca_summary_table"), plotOutput("pca_screeplot"), plotOutput("pca_biplot")),
      tabPanel("Gene Ontology" , verbatimTextOutput("gene_ontology")),
	  tabPanel("Network Analysis", simpleNetworkOutput("networkData", width= "100%", height="1500px"))
    )
  )

)
)