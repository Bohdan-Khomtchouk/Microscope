# Copyright (C) 2015-2016 Bohdan Khomtchouk and James Hennessy

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
library(png)

# backend 
server <- shinyServer(function(input, output) {	
  
  # instructions tab
  output$text1 <- renderText({ "0) You can easily make a .csv file by simply saving your Microsoft Excel workbook as a .csv through 'Save As'.  Before saving as a .csv, your Excel file should look something like:" })
  output$text2 <- renderText({ "Please note that all cell values must be positive (i.e., corresponding to raw gene expression values, i.e., read counts per gene per sample) from a ChIP-seq or RNA-seq experiment.  A sample .csv file is provided under the 'Download Sample Input File' button.  Press that button, save the file to your computer, then click the 'Choose File' button to upload it.  In offbeat cases where the input file is a combination of various standard (or non-standard) delimiters, simply use the 'Text to Columns' feature in Microsoft Excel under the 'Data' tab to parse the file before using MicroScope." })
  output$text3 <- renderText({ "1) After uploading a .csv file, navigate to the DE analysis tab, and specify 'Non-Control Samples' in the box below the button. By default, all remaining samples will be designated as control samples. After pressing the 'Run Statistics' button, your statistical table will appear in the DE Analysis panel. Positive values of log2(FC) (i.e., fold change) indicate upregulation in experimental samples (i.e., experimental samples have higher expression values relative to controls). Negative values indicate downregulation in experimental samples (i.e., control samples have higher expression values relative to experimentals). See the MicroScope publication for more details on the differential expression analysis." })
  output$text4 <- renderText({ "2) To draw a heatmap of these differential expression results, choose whether you want to use either the FDR cutoff or the P-value cutoff or both, by using the selection inputs (this will only impact the heatmap itself).  After clicking the 'Draw Heatmap' button, navigate to the Heatmap panel to see the resultant heatmap. After you see a heatmap, click and drag anywhere in the heatmap to zoom in. Click once to zoom out. Note that this click-drag-zoom action is not recommended for large input datasets (>= 5000 genes), as you may get unexpectedly disconnected from the R Shiny servers performing this computationally intensive JavaScript-based task. In cases like this, it is best to either re-login/re-upload/wait, or pre-filter your input dataset down to a set of genes of specific interest. You may also log-transform your data (and the resultant heatmap) automatically, change color schemes, perform hierarchical clustering, etc.." })
  output$text5 <- renderText({ "3) To perform principal component analysis on your heatmap, specify your matrix type (i.e., covariance or correlation matrix) in the sidebar panel marked 'Choose PCA Option'.  After pressing the 'Run PCA' button, your PCA results will appear in the PCA panel.  You may download both the biplot and screeplot to your computer using the buttons provided.  See the MicroScope publication for more details on the principal component analysis suite." })
  output$text6 <- renderText({ "4) Feel free to download either the heatmap, PCA graphs, and/or differential expression analysis table to your computer using the buttons provided." })  
  output$text7 <- renderText({ "5) You may now navigate to the Gene Ontology panel for detailed instructions about performing GO analysis on the differentially expressed genes.  You may then download the gene ontology results to your computer." })  
  output$text8 <- renderText({ "6) After performing GO analysis, navigate to the Network Analysis panel for information about performing network analysis on the differentially expressed genes.  You may interact with the resultant network by clicking, zooming, and dragging any element of the network.  You may also download the network analysis results to your computer." })  
  output$text9 <- renderText({ "7) For more information about this software, please visit the MicroScope publication.  To run MicroScope on a new dataset, simply refresh the page or open a new browser window." })
  

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
    	need(input$filename != 0, "To perform statistical analysis, please select a file for input") 
    )
    inFile <- input$filename
    if (is.null(inFile)) return(NULL)
    read.table(inFile$datapath, header= TRUE, sep=",", quote='"', row.names=1)
  })
  
  
  # filter stats table based on cutoffs
  slimStats <- reactive({
    results_df <- stats()
    datasetInputTwo <- datasetInput()

    if (input$pvFDRchoose == "Pvalue"){
        pNewResults_df <- results_df[which(results_df$PValue<input$statPV),]
        pNames <- row.names(pNewResults_df)
        pFinalResult <- datasetInputTwo[pNames,]
    }
    else if(input$pvFDRchoose == "FDR"){
        fNewResults_df <- results_df[which(results_df$FDR<input$statFDR),]
        fNames <- row.names(fNewResults_df)
        fFinalResult <- datasetInputTwo[fNames,]
    }
    else if(input$pvFDRchoose == "both"){
        bNewResults_df <- results_df[which(results_df$FDR<input$statFDR & results_df$PValue<input$statPV),]
        bNames <- row.names(bNewResults_df)
        bFinalResult <- datasetInputTwo[bNames,]
    }
  })
  
  
  # heatmap height
    output$pixelation <- renderUI({
     
      slimStats()
      inputLines <- NROW(slimStats())
      if(inputLines >= 0 && inputLines <= 2001){
        d3heatmapOutput("heatmap", width = "100%", height = "700px")
      }
      else if(inputLines > 2001 && inputLines <= 5001){
        d3heatmapOutput("heatmap", width = "100%", height = "1500px")
      }
      else if(inputLines > 5001 && inputLines <= 10001){
        d3heatmapOutput("heatmap", width = "100%", height = "3000px")
      }
      else if(inputLines > 10001 && inputLines <= 20001){
        d3heatmapOutput("heatmap", width = "100%", height = "8000px")
      }
      else if(inputLines > 20001 && inputLines <= 40001){
        d3heatmapOutput("heatmap", width = "100%", height = "12000px")
      }
      else if(inputLines > 40001 && inputLines <= 60001){
        d3heatmapOutput("heatmap", width = "100%", height = "18000px")
      }
      else
        d3heatmapOutput("heatmap", width = "100%", height = "30000px")
      
    })	
  
  
  # log2 data transformation
  log2_datasetInput <- reactive({
    slimStats()
    cpm(slimStats(), prior.count=2, log=TRUE)
  })


  # d3heatmap prep					
  plot <- reactive({
    df = slimStats()
    if (input$goButtonHeat == 0) {return(validate(
      need(input$filename != 0, "To create a heatmap, please first select a file for input"),
      need(input$goButtonHeat !=0 , "To create a heatmap showing statistically significant genes, first specify your parameters and then press 'Draw Heatmap' button"),
      need(!is.null(df), "No heatmap to display at these statistical cutoff parameters")
    ))}   
    else {
    	d3heatmap( 
      		if (input$log2_transformed_data) log2_datasetInput() else slimStats(),
      		cexRow = as.numeric(as.character(input$xfontsize)),
      		cexCol = as.numeric(as.character(input$yfontsize)),
      		colors = input$choose,
      		k_row = input$color_row_branches,
      		k_col = input$color_column_branches,
      		dendrogram = input$dendrogram
    	)  
  	}
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
  
  
  # PCA engine
  PCAfun <- function() {
        df <- datasetInput()
  		rownames(df) <- c()
  		dm <- data.matrix(df)
  		if (input$Type == "Covariance Matrix") {
  			PCA <<- prcomp(dm, scale = TRUE)
  		}  
  		else if (input$Type == "Correlation Matrix") {
  			PCA <<- prcomp(dm, scale = FALSE)
  		}
  	}
  
  
  # PCA biplot
  output$pca_biplot <- renderPlot({
    validate(
    	need(input$pcaButton !=  0, "To generate a biplot, then click 'Run PCA'") 
    	)
    PCAfun()
	biplot(PCA, scale = 0)
	mtext("Biplot", line = 3, col = "black", font = 2, cex = 1.2) 
  })
  
  
  # biplot download								
  output$downloadBiplot <- downloadHandler(
    filename <- function() {
      paste0(basename(file_path_sans_ext(input$filename)), '_biplot', '.png', sep='')
    },
    content <- function(file) {
      if(input$pcaButton != 0) {
      		png(file)
      		tiff(
        		file,
        		width = 2000,
        		height = 2000,
        		units = "px",
        		pointsize = 12,
        		res = 300
      		)
      		biplot(
        		PCA,
        		scale = 0,
        		col = c("red", "blue")
      		)
      		dev.off()
    	}
      else {
        return()
      }
    }
  )
  
  
  # PCA screeplot 
  output$pca_screeplot <- renderPlot({
    validate(
    	need(input$pcaButton != 0, "To generate a screeplot, click 'Run PCA'") 
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
    	screeplot(Screeplot, type = "lines")  #Screeplot auto-title
    }
  })
  
  
  # screeplot download								
  output$downloadScreeplot <- downloadHandler(
    filename <- function() {
      paste0(basename(file_path_sans_ext(input$filename)), '_screeplot', '.png', sep='')
    },
    content <- function(file) {
      if(input$pcaButton != 0) {
      	png(file)
      	tiff(
        	file,
        	width = 4000,
        	height = 2000,
        	units = "px",
        	pointsize = 12,
        	res = 300
      		)
		screeplot(
			Screeplot, 
			type = "lines"
			)
      	dev.off()
      }
      else {
        return()
      }
    }
  )
  
  
  # PCA summary info
  output$pca_summary_table <- renderPrint({
    if (input$pcaButton == 0) {return(validate(
      need(input$filename != 0, "Table of PCA summary: To conduct principal component analysis, please first select a file for input") %then%
        need(input$goButton != 0, "Table of PCA summary: To conduct principal component analysis, click 'Run PCA'")
    ))}        
    else {
		PCAfun()
    	summary(PCA)
    }
  }) 
  
  
  # edgeR prep
  output$expcolumns <- renderUI({
    df <- datasetInput()
    if (is.null(df)) return(NULL)
    expcolumns <- names(df)
    selectInput('expcolumns', 'Specify Non-Control Samples:', choices = expcolumns, multiple = TRUE)
  })			
  
  
  # edgeR statistical engine
  stats <- reactive({
    if(!is.null(datasetInput()) & !is.null(input$expcolumns)){
   
    group <- as.numeric(names(datasetInput()) %in% input$expcolumns)
    y <- DGEList(counts = datasetInput(), group = group)
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    et <- exactTest(y)
    results <- topTags(et, n=50000)
    results_df <- as.data.frame(results)
      
    }
  })
  
  
  # differential expression analysis table
  `%then%` <- shiny:::`%OR%`
  output$table <- renderTable({
    if (input$goButton == 0) {return(validate(
    	need(input$filename != 0, "To run differential expression analysis, please first select a file for input") %then%
        need(input$goButton != 0, "To run differential expression analysis, please select your experimental samples under 'Specify Non-Control Samples' then click 'Run Statistics'")
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

    group <<- as.numeric(names(datasetInput()) %in% input$expcolumns)
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
        require("org.Mm.eg.db")
      }
      else if(input$Genome == "hg19"){
        require("org.Hs.eg.db")
      }
      else if(input$Genome == "danRer6"){
        require("org.Dr.eg.db")
      }
      else if(input$Genome == "ce6"){
        require("org.Ce.eg.db")  
      }
      else if(input$Genome == "panTro2"){
        require("org.Pt.eg.db")  
      }
      else if(input$Genome == "rn4"){
        require("org.Rn.eg.db")  
      }
      else if(input$Genome == "dm3"){
        require("org.Dm.eg.db")  
      }
      else if(input$Genome == "sacCer2"){
        require("org.Sc.sgd.db")  
      }
      else if(input$Genome == "bosTau4"){
        require("org.Bt.eg.db")  
      }
      else if(input$Genome == "canFam2"){
        require("org.Cf.eg.db")  
      }
      else if(input$Genome == "anoGam1"){
        require("org.Ag.eg.db")  
      }
      else if(input$Genome == "rheMac2"){
        require("org.Mmu.eg.db")  
      }
      else if(input$Genome == "xenTro2"){
        require("org.Xl.eg.db")  
      }
      else if(input$Genome == "galGal3"){
        require("org.Gg.eg.db")  
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
      need(input$goButton != 0, "For gene ontology analysis, please first select your experimental samples under 'Specify Non-Control Samples' then click 'Run Statistics'") %then%
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
      need(input$goButton != 0, "For network analysis, please first select your experimental samples under 'Specify Non-Control Samples' then click 'Run Statistics'") %then%
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