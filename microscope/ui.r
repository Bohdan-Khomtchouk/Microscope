# Copyright (C) 2015-2016 Bohdan Khomtchouk, Vytas Dargis-Robinson, and James Hennessy

# This file is part of MicroScope.

# MicroScope is an R Shiny and JavaScript (D3.js) software program designed to produce dynamic, interactive heatmaps in an Internet browser.
# MicroScope allows you to magnify any portion of a heatmap by a simple click-and-drag feature to zoom in, and click-once feature to zoom out.
# MicroScope is designed with large heatmaps in mind (e.g., gene expression heatmaps with thousands of genes), where individual entries quickly become unreadable as more are added. 
# However, MicroScope allows you to repeatedly zoom in to any sector of the heatmap to investigate a region, cluster, or even a single gene.  
# MicroScope also allows you to hover the mouse pointer over any specific gene to show precise expression level details.
# In addition to visual magnification, MicroScope also allows users to analytically perform real-time statistical analyses, as well as perform gene ontology and network analyses, with simple button clicks.

# For more information, please see: "Khomtchouk BB, Dargis-Robinson V, Hennessy JR, Wahlestedt C. "MicroScope: real-time statistical analytics and visualization platform for gene expression heatmaps". bioRxiv doi: http://dx.doi.org/10.1101/034694"

# MicroScope is an ongoing bioinformatics software project fully financially supported by the United States Department of Defense (DoD) 
# through the National Defense Science and Engineering Graduate Fellowship (NDSEG) Program. This research was conducted with Government support 
# under and awarded by DoD, Army Research Office (ARO), National Defense Science and Engineering Graduate (NDSEG) Fellowship, 32 CFR 168a.

# Please cite: "Khomtchouk BB, Dargis-Robinson V, Hennessy JR, Wahlestedt C. "MicroScope: real-time statistical analytics and visualization platform for gene expression heatmaps". bioRxiv doi: http://dx.doi.org/10.1101/034694" within any source that makes use of any methods inspired by MicroScope.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,c
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
    sliderInput("buffer", "Buffer Size:", value = 700, min = 700, max = 20000),
    checkboxInput("log2_transformed_data", "log2 Transform Heatmap"),
    selectInput("choose", "Choose Color Scheme:", c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu", "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges", "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues")),
    selectInput("dendrogram", "Apply Clustering:", c("none", "row", "column", "both")),
    numericInput("color_row_branches", "Color Row Branches:", value = 1),
    numericInput("color_column_branches", "Color Column Branches:", value = 1),
    sliderInput("xfontsize", "Choose X Font Size:", min = 0.3, max = 2, value = 0.5),
    sliderInput("yfontsize", "Choose Y Font Size:", min = 0.3, max = 2, value = 1.0),
    downloadButton("downloadHeatmap", "Download Heatmap"),
    uiOutput("ctrlcolumns"),
    actionButton("goButton", "Run Statistics"),
    downloadButton("downloadtable", "Download Stats Table")
    ,
    conditionalPanel(
      condition = "input.goButton == true",
      selectInput("Genome", "Choose Genome Database:", selected = "mm9", c("Mouse" = "mm9", "Human" = "hg38", "Chimpanzee" = "panTro4", "Rat" = "rn6", "Worm" = "ce10", "Zebrafish" = "danRey10", "Fly" = "dm6", "Yeast" = "sacCer3", "Frog" = "xenTro7", "Chicken" = "galGal4", "Pig" = "susScr3", "Cow" = "bosTau8", "Dog" = "canFam3")),
      selectInput("geneRef", "Choose Gene Identifier:", selected = "geneSymbol", c("Gene Symbol" = "geneSymbol", "Ensembl ID" = "ensGene")),
      numericInput("numberGenes", "Choose How Many Top Gene Ontologies to Display:", value = 10),
      selectInput("chooseEnriched", "Stratify Top Gene Ontologies By:", c("P-value", "FDR")),
      sliderInput("cutoffP", "Choose Gene Ontology P-value Cutoff:", min = 0.001, max = 0.1, value = 0.05),
	  sliderInput("cutoffFdr", "Choose Gene Ontology FDR Cutoff:", min = 0.001, max = 0.1, value = 0.05),
	  actionButton("goData","Do Gene Ontology Analysis")
    ),
    conditionalPanel(
      condition = "input.goData == true",
      actionButton("doNets", "Do Network Analysis"),
      downloadButton("downloadSimpleNetwork", "Download Network")
    )
  ),


  mainPanel(
    tabsetPanel(
      tabPanel("Instructions", textOutput("text1"), img(src='excel.png'), textOutput("text2"), textOutput("text3"), textOutput("text4"), textOutput("text5"), textOutput("text6"), textOutput("text7"), textOutput("text8")),
      tabPanel("Heatmap", uiOutput("heatmapOutput")),
      tabPanel("Statistical Analysis", tableOutput("table")),
      tabPanel("Gene Ontology" , verbatimTextOutput("gene_ontology")),
	  tabPanel("Network Analysis", simpleNetworkOutput("networkData", width= "100%", height="1500px"))
    )
  )

)
)