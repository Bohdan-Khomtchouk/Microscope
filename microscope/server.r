# Copyright (C) 2015 Bohdan Khomtchouk, James Hennessy, and Vytas Dargis-Robinson

# This file is part of MicroScope.

# MicroScope is an RShiny and JavaScript (D3.js) software program designed to produce dynamic, interactive heatmaps in an Internet browser.
# MicroScope allows you to magnify any portion of a heatmap by a simple click-and-drag feature to zoom in, and click-once feature to zoom out.
# MicroScope is designed with large heatmaps in mind (e.g., gene expression heatmaps with thousands of genes), where individual entries quickly become unreadable as more are added. 
# However, MicroScope allows you to repeatedly zoom in to any sector of the heatmap to investigate a region, cluster, or even a single gene.  
# MicroScope also allows you to hover the mouse pointer over any specific gene to show precise expression level details.

# MicroScope is an ongoing bioinformatics software project fully financially supported by the United States Department of Defense (DoD) 
# through the National Defense Science and Engineering Graduate Fellowship (NDSEG) Program. This research was conducted with Government support 
# under and awarded by DoD, Army Research Office (ARO), National Defense Science and Engineering Graduate (NDSEG) Fellowship, 32 CFR 168a.

# Current work is underway to expand Microscope's user-friendly features.

# Please cite: "Khomtchouk et al.: 'MicroScope: magnifying interactive heatmaps using RShiny and JavaScript', 2015 (in preparation)" 
# within any source that makes use of any methods inspired by MicroScope. 

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

library(d3heatmap)
library(shiny)



shinyServer(function(input, output) {
  output$contents <- renderTable({
    pwd<-"/path/to/my_directory/MICROSCOPE/genes_file.csv"
genes.numeric<-read.csv(pwd, header= TRUE, sep=",", quote= '"',row.names=1)
    # input$file1 will be NULL initially. After the user selects and uploads a 
    # file, it will be a data frame with 'name', 'size', 'type', and 'datapath' 
    # columns. The 'datapath' column will contain the local filenames where the 
    # data can be found.

data <- reactive({
    dist <- switch(input$dist,
                   norm = rnorm,
                   unif = runif,
                   lnorm = rlnorm,
                   exp = rexp,
                   rnorm)
    
    dist(input$n)
  })
   output$summary <- renderPrint({
     dataset <- datasetInput()

    capture.output(summary(dataset),file="genes.numeric")

  })
  
   	output$bake<-renderPrint({(genes.numeric[row.names(genes.numeric)=="geneA ",])})

	datasetInput<- reactive({
   	output$something<-renderPrint({(genes.numeric[row.names(genes.numeric)=="geneA ",])})
      output$names <- input$text
	output$something2<-renderPlot({plot(genes.numeric$entry1,genes.numeric$entry2)})
	output$color<-input$color
 	inFile <- input$file1
    	if (is.null(inFile))return(NULL)
     percent<-input$percentage
	output$table<-read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)})
    output$rtable<-renderDataTable( table())
 
  })
})