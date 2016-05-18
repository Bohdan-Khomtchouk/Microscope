# Copyright (C) 2015-2016 Bohdan Khomtchouk, Vytas Dargis-Robinson, and James Hennessy

# This file is part of MicroScope.

# MicroScope is an R Shiny and JavaScript (D3.js) software program designed to produce dynamic, interactive heatmaps in an Internet browser.
# MicroScope allows you to magnify any portion of a heatmap by a simple click-and-drag feature to zoom in, and click-once feature to zoom out.
# MicroScope is designed with large heatmaps in mind (e.g., gene expression heatmaps with thousands of genes), where individual entries quickly become unreadable as more are added. 
# However, MicroScope allows you to repeatedly zoom in to any sector of the heatmap to investigate a region, cluster, or even a single gene.  
# MicroScope also allows you to hover the mouse pointer over any specific gene to show precise expression level details.
# In addition to visual magnification, MicroScope also allows users to analytically perform real-time statistical analyses, as well as perform gene ontology and network analyses, with simple button clicks.

# For more information, please see: "Khomtchouk BB, Dargis-Robinson V, Hennessy JR, Wahlestedt C. "MicroScope: real-time statistical analytics and visualization platform for gene expression heatmaps". bioRxiv doi: http://dx.doi.org/10.1101/034694"

# MicroScope is an ongoing bioinformatics software project financially supported by the United States Department of Defense (DoD) 
# through the National Defense Science and Engineering Graduate Fellowship (NDSEG) Program. This research was conducted with Government support 
# under and awarded by DoD, Army Research Office (ARO), National Defense Science and Engineering Graduate (NDSEG) Fellowship, 32 CFR 168a.

# Please cite: "Khomtchouk BB, Dargis-Robinson V, Hennessy JR, Wahlestedt C. "MicroScope: real-time statistical analytics and visualization platform for gene expression heatmaps". bioRxiv doi: http://dx.doi.org/10.1101/034694" within any source that makes use of any methods inspired by MicroScope.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ------------------------------------------------------------------------------------

library(shiny)
library(d3heatmap)
library(htmlwidgets)
library(tools)
library(edgeR)
library(GO.db)
library(goseq)
library(networkD3)
library(data.table)
library(dplyr)


# backend 
server <- shinyServer(function(input, output) {	
  
  # instructions tab
  output$text1 <- renderText({ "0) You can easily make a .csv file by simply saving your Microsoft Excel workbook as a .csv through 'Save As'.  Before saving as a .csv, your Excel file should look something like:" })
  output$text2 <- renderText({ "Please note that all cell values must be positive (i.e., corresponding to raw gene expression values, i.e., read counts per gene per sample).  After you upload your .csv file, you can log-transform your data to normalize it (and the resultant heatmap) automatically within MicroScope." })
  output$text3 <- renderText({ "1) After uploading a .csv file, navigate to the Heatmap panel to see the resultant heatmap.  Click and drag anywhere in the heatmap to zoom in.  Click once to zoom out." })
  output$text4 <- renderText({ "2) To perform statistical analysis on your heatmap, specify your control samples in the sidebar panel marked 'Specify Control Samples:' (by default, the remaining samples will be designated as experimental samples).  After pressing the 'Run Statistics' button, your statistical table will appear in the Statistical Analysis panel. Positive values of log(FC) (i.e., fold change) indicate upregulation in experimental samples (i.e., experimental samples have higher expression values relative to controls). Negative values indicate downregulation in experimental samples (i.e., control samples have higher expression values relative to experimentals). " })
  output$text5 <- renderText({ "3) Feel free to download either the heatmap or the statistical analysis table to your computer using the buttons provided." })  
  output$text6 <- renderText({ "4) After generating a heatmap and running statistics on it, navigate to the Gene Ontology panel for information about performing GO analysis on the top differentially expressed genes in your heatmap.  You may then download the gene ontology results to your computer." })  
  output$text7 <- renderText({ "5) After performing GO analysis, navigate to the Network Analysis panel for information about performing network analysis on the top differentially expressed genes in your heatmap.  You may interact with the resultant network by clicking, zooming, and dragging any element of the network.  You may also download the network analysis results to your computer." })  
  output$text8 <- renderText({ "6) For more information about this software, please visit the MicroScope publication." })
  
  
  # sample file download
  output$downloadData <- downloadHandler(
  	filename <- function() {
    	paste('genes', '_file', '.csv', sep='')
  	},
  	content <- function(file) {
    	file.copy("genes_file.csv", file)
  	},
  	contentType = "text/csv"
	)
	
  
  # file upload
  datasetInput <- reactive({
    validate(
    	need(input$filename != 0, "To generate a heatmap, please select a file for input") 
    )
    inFile <- input$filename
    if (is.null(inFile)) return(NULL)
    read.table(inFile$datapath, header= TRUE, sep=",", quote='"', row.names=1)
  })	
  
  
  # log2 data transformation
  log2_datasetInput <- reactive({
    cpm(datasetInput(), prior.count=2, log=TRUE)
  })


  # d3heatmap prep					
  plot <- reactive({
    d3heatmap( 
      if (input$log2_transformed_data) log2_datasetInput() else datasetInput(),
      cexRow = as.numeric(as.character(input$xfontsize)),
      cexCol = as.numeric(as.character(input$yfontsize)),
      colors = input$choose,
      k_row = input$color_row_branches,
      k_col = input$color_column_branches,
      dendrogram = input$dendrogram
    )  
  })
  
  
  # d3heatmap output
  output$heatmap <- renderD3heatmap({
    if(!is.null(datasetInput()))
      plot()
  })
  
  
  # d3heatmap download								
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste0(basename(file_path_sans_ext(input$filename)), '.html')
    },
    content = function(file) {
      saveWidget(plot(), file)
    }
  )
  
  
  # edgeR prep
  output$ctrlcolumns <- renderUI({
    df <- datasetInput()
    if (is.null(df)) return(NULL)
    ctrlcolumns <- names(df)
    selectInput('ctrlcolumns', 'Specify Control Samples:', choices = ctrlcolumns, multiple = TRUE)
  })			
  
  
  # edgeR statistical engine
  stats <- reactive({
    group <- as.numeric(names(datasetInput()) %in% input$ctrlcolumns)
    y <- DGEList(counts = datasetInput(), group = group)
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    et <- exactTest(y)
    results <- topTags(et, n=50000)
    results_df <- as.data.frame(results)
  })


  # statistical analysis table
  `%then%` <- shiny:::`%OR%`
  output$table <- renderTable({
    if (input$goButton == 0) {return(validate(
    	need(input$filename != 0, "To run statistics, please first select a file for input") %then%
        need(input$goButton != 0, "To run statistics, please select your control samples under 'Specify Control Samples' then click 'Run Statistics'")
    ))}        
    else {
      stats()
    }
  }, digits = -3)
  
  
  # statistical table download
  output$downloadtable <- downloadHandler(
    filename = function() {
      paste(basename(file_path_sans_ext(input$filename)), '_edgeR_stats_table', '.csv', sep='')
    },
    content = function(file) {
      if (input$goButton != 0) write.csv(stats(), file) else return()
    }
  )

   
  # create enriched gene database from GO analysis
  enriched <- reactive({

    group <<- as.numeric(names(datasetInput()) %in% input$ctrlcolumns)
    y <<- DGEList(counts = datasetInput(), group = group)
    y <<- calcNormFactors(y)
    y <<- estimateCommonDisp(y)
    y <<- estimateTagwiseDisp(y)
    et <<- exactTest(y)
    
    genes <<- as.integer(p.adjust(et$table$PValue[et$table$logFC!=0], method="BH") < input$cutoffP)
    genesP <<- as.integer(p.adjust(et$table$PValue, method="BH") < input$cutoffP)
    names(genes) = row.names(et$table[et$table$logFC!=0,])
    pwf <<- nullp(genes, input$Genome, input$geneRef)
    pwf_1 <<- pwf[which(pwf[,1] == 1),]
    rownamespwfOne <<- row.names(pwf_1)
    
	callData <<- reactive({ 
      if(input$Genome == "mm9"){
        library("org.Mm.eg.db")
      }
      else if(input$Genome == "hg38"){
        library("org.Hs.eg.db")
      }
      else if(input$Genome == "danRey10"){
        library("org.Dr.eg.db")
      }
      else if(input$Genome == "ce10"){
        library("org.Ce.eg.db")  
      }
      else if(input$Genome == "panTro4"){
        library("org.Pt.eg.db")  
      }
      else if(input$Genome == "rn6"){
        library("org.Rn.eg.db")  
      }
      else if(input$Genome == "dm6"){
        library("org.Dm.eg.db")  
      }
      else if(input$Genome == "sacCer3"){
        library("org.Sc.sgd.db")  
      }
      else if(input$Genome == "xenTro7"){
        library("MeSH.Xtr.eg.db")  
      }
      else if(input$Genome == "galGal4"){
        library("MeSH.Gga.eg.db")  
      }
      else if(input$Genome == "susScr3"){
        library("org.Ss.eg.db")  
      }
      else if(input$Genome == "bosTau8"){
        library("org.Bt.eg.db")  
      }
      else if(input$Genome == "canFam3"){
        library("org.Cf.eg.db")  
      }
    })
  
    GO.wall <<- goseq(pwf, input$Genome, input$geneRef)
    
    if(input$chooseEnriched == "FDR") ({
    	FDR.enriched.GO <<- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH") < input$cutoffFdr]
    })
    else if(input$chooseEnriched == "P-value") ({
    	enriched.GO <<- GO.wall$category[(GO.wall$over_represented_pvalue) < input$cutoffP]
    })
    
  })
  
  
  # output the enriched genes
  output$gene_ontology <- renderPrint({
  
    if (input$goData == 0) {return(validate(
      need(input$filename != 0, "For gene ontology analysis, please first select a file for input") %then%
      need(input$goButton != 0, "For gene ontology analysis, please first select your control samples under 'Specify Control Samples' then click 'Run Statistics'") %then%
      need(input$goData != 0, "For gene ontology analysis, please select a genome under 'Choose Genome Database', choose the appropriate gene identifier, choose how many top gene ontologies to display, choose how to stratify your top gene ontologies (e.g., by p-value or by FDR), choose the respective statistical cutoff, and then click 'Do Gene Ontology Analysis'.  Now wait while your genome-wide annotations are automatically downloaded from UCSC (this may take, on average, between 15 seconds and 2 minutes, depending on your network connection).  Your gene ontology results will appear shortly.")
    ))}        
    else {
      enriched()
      if(input$chooseEnriched == "FDR") ({
        enriched.GO = FDR.enriched.GO   
      })
    }
    
  if(!length(enriched.GO) == 0) {
    for(go in enriched.GO[1:input$numberGenes]) {
      print(GOTERM[[go]])
      cat("--------------------------------------\n")
    }
  }
  else{
    print("No enriched genes present at this statistical cutoff level")
  }
  
  })
  
  
  # simple network engine prep
  grow <- function(nominal_GO_df, rownamespwfOne, enriched.GO, Genome, topEnriched) ({
    total = list()
    getGo <- getgo(rownamespwfOne, Genome, input$geneRef)
    Base = data.table(rownamespwfOne, getGo)
    for(go in nominal_GO_df$ind[1:topEnriched]){
        Names <- Base$rownamespwfOne[grep(go,Base$getGo)]
        nameList = list(Names) 
        total = append(total, nameList)
    }
    names(total) <- nominal_GO_df$ind[1:topEnriched]
    total <- total
  })
  
  
  # simple network engine
  network <- reactive({
    
    if(input$chooseEnriched == "FDR")({
    	enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH") < input$cutoffFdr]
    })
    else if(input$chooseEnriched == "P-value")({
    	enriched.GO <- GO.wall$category[(GO.wall$over_represented_pvalue) < input$cutoffP]
    })
    
    
    #network analysis loop (FDR vs nominal p-value) within simple network engine
    if(!length(enriched.GO) == 0){
    	nominal_GO_df <- stack(lapply(mget(enriched.GO, GOTERM, ifnotfound = NA), Term)) 
    
    	total <- grow(nominal_GO_df, rownamespwfOne, input$chooseEnriched, input$Genome, input$numberGenes)
    	namesGo <- names(total)
    	vector = c()
    	targetVector = c()
    
    	for (go in 1:length(namesGo)){
    		lengthName <- length(total[[go]])
    		for(populate in 1:lengthName){
    			nameGene <- nominal_GO_df$values[grep(namesGo[go], nominal_GO_df$ind)]
      			vector = append(vector, nameGene)
      			targetVector = append(targetVector, total[[go]][populate])
      		}
    	}
    
    	src <- vector
    	target <- targetVector
    	networkData <- data.frame(src, target)   
    	simpleNetwork(networkData, charge = -1000, opacity = 1, zoom = T)
    }
    
	else {
		return(validate(need(input$goData == 0, "No enriched genes are available to make a network at this statistical cutoff level")))
    }
    
  })
  
  
  # simple network output
  output$networkData <- renderSimpleNetwork({
	if (input$doNets == 0) {return(validate(
      need(input$filename != 0, "For network analysis, please first select a file for input") %then%
      need(input$goButton != 0, "For network analysis, please first select your control samples under 'Specify Control Samples' then click 'Run Statistics'") %then%
      need(input$goData != 0, "For network analysis, please first select a genome under 'Choose Genome Database' then click 'Do Gene Ontology Analysis'") %then%
      need(input$doNets != 0, "For network analysis, please click on 'Do Network Analysis' and wait for about 5-30 seconds")
    ))}
	else {
      network()
    }
  })
  
  
  # simple network download								
  output$downloadSimpleNetwork <- downloadHandler(
    filename = function() {
      paste0(basename(file_path_sans_ext(input$filename)), '_network', '.html', sep='')
    },
    content = function(file) {
    	if (input$doNets != 0) {
    		saveWidget(network(), file)
    	} 
    	else {
    		return()
    	}
    }
  )
  

  
})