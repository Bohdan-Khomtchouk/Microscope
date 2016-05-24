# Copyright (C) 2015-2016 Bohdan Khomtchouk, James Hennessy, and Vytas Dargis-Robinson

# This file is part of MicroScope.

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
  output$text2 <- renderText({ "Please note that all cell values must be positive (i.e., corresponding to raw gene expression values, i.e., read counts per gene per sample) from a ChIP-seq or RNA-seq experiment.  Users with microarray data are advised to only use MicroScope’s heatmap utility (but not the statistical/GO/network analysis utilities).  A sample .csv file is provided under the 'Download Sample Input File' button.  Press that button, save the file to your computer, then click the 'Choose File' button to upload it.  In offbeat cases where the input file is a combination of various standard (or non-standard) delimiters, simply use the 'Text to Columns' feature in Microsoft Excel under the 'Data' tab to parse the file before using MicroScope." })
  output$text3 <- renderText({ "1) After uploading a .csv file, navigate to the Heatmap panel to see the resultant heatmap.  If your dataset is extremely big and you do not see a heatmap (or it appears empty), simply drag the 'Buffer Size' button and MicroScope will adapt and allocate resources to your big dataset.  After you see a heatmap, click and drag anywhere in the heatmap to zoom in.  Click once to zoom out.  You may also log-transform your data (and the resultant heatmap) automatically within MicroScope." })
  output$text4 <- renderText({ "2) To perform statistical analysis on your heatmap, specify your control samples in the sidebar panel marked 'Specify Control Samples'.  By default, all remaining samples will be designated as experimental samples.  After pressing the 'Run Statistics' button, your statistical table will appear in the Statistical Analysis panel. Positive values of log2(FC) (i.e., fold change) indicate upregulation in experimental samples (i.e., experimental samples have higher expression values relative to controls). Negative values indicate downregulation in experimental samples (i.e., control samples have higher expression values relative to experimentals).  See the MicroScope publication for more details on the statistical analysis. " })
  output$text5 <- renderText({ "3) Feel free to download either the heatmap or the statistical analysis table to your computer using the buttons provided." })  
  output$text6 <- renderText({ "4) After generating a heatmap and running statistics on it, navigate to the Gene Ontology panel for detailed instructions about performing GO analysis on the top differentially expressed genes in your heatmap.  You may then download the gene ontology results to your computer." })  
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

  
  # set buffer size for big data users
  output$heatmapOutput <- renderUI({
	d3heatmapOutput("heatmap", height = paste0(input$buffer, "px"))
  })


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
  
  
  # Unlock features (PCA and stats UI) after heatmap upload
  outputOptions(output, 'heatmap', suspendWhenHidden = FALSE)
  
  
  # d3heatmap download								
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste0(basename(file_path_sans_ext(input$filename)), '.html')
    },
    content = function(file) {
      saveWidget(plot(), file)
    }
  )
  
  
  # PCA biplot
  output$pca_biplot <- renderPlot({
    validate(
    	need(input$pcaButton !=  0, "To generate a biplot, follow the instructions above, then click 'Run PCA'") 
    	)
    if (length(PCA) != 0){
    	df <- datasetInput()
  		rownames(df) <- c()
  		dm <- data.matrix(df)
  		if (input$Type == "Covariance Matrix") {
  			PCA <<- prcomp(dm, scale = TRUE)
  		}  
  		else if (input$Type == "Correlation Matrix") {
  			PCA <<- prcomp(dm, scale = FALSE)
  		}
		biplot(PCA, scale = 0)
		mtext("Biplot", line = 3, col = "black", font = 2, cex = 1.2)
    } 
  })
  
  
  # PCA screeplot 
  output$pca_screeplot <- renderPlot({
    validate(
    	need(input$pcaButton != 0, "To generate a screeplot, follow the instructions above, then click 'Run PCA'") 
    	)
    if(length(PCA) != 0){
    	df <- datasetInput()
  		rownames(df) <- c()
  		dm <- data.matrix(df)
  		if (input$Type == "Covariance Matrix") {
  			Screeplot <<- prcomp(dm, scale = TRUE)
  		}  
  		else if (input$Type == "Correlation Matrix") {
  			Screeplot <<- prcomp(dm, scale = FALSE)
  		}
    	screeplot(Screeplot, type = "lines")
    }
  })
  
  
  # PCA summary info
  output$pca_summary_table <- renderPrint({
    if (input$pcaButton == 0) {return(validate(
      need(input$filename != 0, "Table of PCA summary: To conduct principal components analysis, please first select a file for input") %then%
        need(input$goButton != 0, "Table of PCA summary: To conduct principal components analysis, click 'Run PCA'")
    ))}        
    else {
    	df <- datasetInput()
  		rownames(df) <- c()
  		dm <- data.matrix(df)
  		if (input$Type == "Covariance Matrix") {
  			PCA <<- prcomp(dm, scale = TRUE)
  		}  
  		else if (input$Type == "Correlation Matrix") {
  			PCA <<- prcomp(dm, scale = FALSE)
  		}
    	summary(PCA)
    }
  }) 
  
  
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
