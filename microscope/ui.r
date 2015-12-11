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
library(RColorBrewer)

pwd<-"C:/Users/James.JAMES/Desktop/Microscope-master/genes_file.csv"
genes<-read.csv(pwd, header= TRUE, sep=",", quote= '"',row.names=1)
genes.numeric <- data.matrix(genes)


shinyUI(pageWithSidebar(
  headerPanel("MicroScope"),
  sidebarPanel(
    fileInput('file1', 'Choose CSV File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    tags$hr(),
    checkboxInput('header', 'Header', TRUE),
    radioButtons('sep', 'Separator',
                 c(Comma=',',
                   Semicolon=';',
                   Tab='\t'),
                 ','),

    br(),
    

    radioButtons('quote', 'Quote',
                 c(None='',
                   'Double Quote'='"',
                   'Single Quote'="'"),
                 '"'),
    radioButtons('color', 'Colors',
                 c(Blue ='Blues',
                   Greens='Greens'),
                 'Blues'),
    radioButtons('percentage', 'Percentage',
                 c(None='',
                   topOne ='.01',
                   topFive='.05',
                   topTen='.1'),
                 'none'),
    textInput("text", label=h3("Text Input"), value="enter Speciman...")

  ),
    
  
  
  mainPanel(
	tabsetPanel(
	
	type="tabs",
	
	tabPanel("Interactive! Heatmap", d3heatmap(genes.numeric, dendrogram ="row", color= colorRamp(c("green", "black", "red"), space = "rgb"), width = "150%", height = "1000px", cexRow=0.5)),
	
	tabPanel("Instructions",textInput("text", label=h3("Click-and-drag to zoom in. Click-once to zoom out."), value="It's that easy!"))
)
)
  ))
