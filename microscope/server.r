# Copyright (C) 2015-2016 Bohdan Khomtchouk, Vytas Dargis-Robinson, and James Hennessy

# This file is part of MicroScope.

# MicroScope is an RShiny and JavaScript (D3.js) software program designed to produce dynamic, interactive heatmaps in an Internet browser.
# MicroScope allows you to magnify any portion of a heatmap by a simple click-and-drag feature to zoom in, and click-once feature to zoom out.
# MicroScope is designed with large heatmaps in mind (e.g., gene expression heatmaps with thousands of genes), where individual entries quickly become unreadable as more are added. 
# However, MicroScope allows you to repeatedly zoom in to any sector of the heatmap to investigate a region, cluster, or even a single gene.  
# MicroScope also allows you to hover the mouse pointer over any specific gene to show precise expression level details.
# In addition to visual magnification, MicroScope also allows users to analytically perform real-time statistical analyses with simple button clicks.

# For more information, please see: "Khomtchouk BB, Dargis-Robinson V, Hennessy JR, Wahlestedt C. “MicroScope: real-time statistical analytics and visualization platform for gene expression heatmaps”. bioRxiv doi: http://dx.doi.org/10.1101/034694"

# MicroScope is an ongoing bioinformatics software project fully financially supported by the United States Department of Defense (DoD) 
# through the National Defense Science and Engineering Graduate Fellowship (NDSEG) Program. This research was conducted with Government support 
# under and awarded by DoD, Army Research Office (ARO), National Defense Science and Engineering Graduate (NDSEG) Fellowship, 32 CFR 168a.

# Please cite: "Khomtchouk BB, Dargis-Robinson V, Hennessy JR, Wahlestedt C. “MicroScope: real-time statistical analytics and visualization platform for gene expression heatmaps”. bioRxiv doi: http://dx.doi.org/10.1101/034694" within any source that makes use of any methods inspired by MicroScope.

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


# backend 
server <- shinyServer(function(input, output) {	

    #instructions tab
    output$text1 <- renderText({ "0) You can easily make a .csv file by simply saving your Microsoft Excel workbook as a .csv through 'Save As'.  Before saving as a .csv, your Excel file should look something like:" })
    output$text2 <- renderText({ "1) Upload a .csv file.  Navigate to the Heatmap panel to see the resultant heatmap.  Click and drag anywhere in the heatmap to zoom in.  Click once to zoom out." })
    output$text3 <- renderText({ "2) To perform statistical analysis on your heatmap AFTER you've loaded in the .csv file, specify your control samples in the sidebar panel marked 'Specify Control Samples:' (by default, the remaining samples will be designated as experimental samples).  After pressing the Do Stats! button, your statistical table will appear in the Statistical Analysis panel." })
    output$text4 <- renderText({ "3) Feel free to download either the heatmap or the statistical analysis table to your computer using the buttons provided.  For more information, please visit the MicroScope publication." })


	#file upload
	datasetInput <- reactive({
    	inFile <- input$filename
    	if (is.null(inFile)) return(NULL)
    	read.table(inFile$datapath, header= TRUE, sep=",", quote='"', row.names=1)
							})	
							
	
	#edgeR prep
	output$ctrlcolumns <- renderUI({
    	df <- datasetInput()
    	if (is.null(df)) return(NULL)
    	ctrlcolumns <- names(df)
    	selectInput('ctrlcolumns', 'Specify Control Samples:', choices = ctrlcolumns, multiple = TRUE)
    })			
	
	
	#edgeR statistical engine
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
	
	
	#statistical analysis table
	output$table <- renderTable({
		if (input$goButton == 0) {return()}        
    	else{
        	stats()
        	}
		}, digits = -3)
	
	
	#statistical table download
	output$downloadtable <- downloadHandler(
    	filename = function() {
    		paste(basename(file_path_sans_ext(input$filename)), '_edgeR_stats_table', '.csv', sep='')
    	},
    	content = function(file) {
    		if (input$goButton != 0) write.csv(stats(), file) else return()
    		}
    	)

	
	#d3heatmap prep					
	plot <- reactive({
 		d3heatmap( 
      		datasetInput(),
      		cexRow = as.numeric(as.character(input$xfontsize)),
      		cexCol = as.numeric(as.character(input$yfontsize)),
      		colors = input$choose,
      		k_row = input$color_row_branches,
      		k_col = input$color_column_branches,
      		dendrogram = input$dendrogram
    			 )  
    })
    
	
	#d3heatmap output
	output$heatmap <- renderD3heatmap({
		if(!is.null(datasetInput()))
        	plot()
    })
    
  	
  	#d3heatmap download								
    output$downloadHeatmap <- downloadHandler(
        filename = function() {
            paste0(basename(file_path_sans_ext(input$filename)), '.html')
        },
        content = function(file) {
            saveWidget(plot(), file)
        }
    )
    
})