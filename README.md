# MicroScope

## About

MicroScope is an R Shiny and JavaScript (D3.js) software program designed to produce dynamic, interactive heatmaps in an Internet browser.  MicroScope allows you to magnify any portion of a heatmap by a simple click-and-drag feature to zoom in, and click-once feature to zoom out.  MicroScope is designed with large heatmaps in mind (e.g., gene expression heatmaps with thousands of genes), where individual entries quickly become unreadable as more are added.  However, MicroScope allows you to repeatedly zoom in to any sector of the heatmap to investigate a region, cluster, or even a single gene.  MicroScope also allows you to hover the mouse pointer over any specific gene to show precise expression level details.  In addition to visual magnification, MicroScope also allows users to analytically perform real-time statistical analyses, as well as perform gene ontology and network analyses, with simple button clicks.

MicroScope is an ongoing bioinformatics software project financially supported by the United States Department of Defense (DoD) through the National Defense Science and Engineering Graduate Fellowship (NDSEG) Program. This research was conducted with Government support under and awarded by DoD, Army Research Office (ARO), National Defense Science and Engineering Graduate (NDSEG) Fellowship, 32 CFR 168a.

Please cite: "Khomtchouk BB, Dargis-Robinson V, Hennessy JR, Wahlestedt C.  “MicroScope: comprehensive genome analysis software suite for gene expression heatmaps” (submitted). bioRxiv doi: http://dx.doi.org/10.1101/034694" within any source that makes use of any methods inspired by MicroScope. 

## Usage (for general public)

##### Just click here!... https://microscope.shinyapps.io/microscope

## Installation (for developers only)

### Requirements for developers

* R programming language
  * RStudio

## How to run (for developers only)

##### Git clone this repo to your computer, and in RStudio type:
* `setwd("~/path/to/my_directory/Microscope")`
* `install.packages("shiny")`
* `library(shiny)`
* `install.packages("d3heatmap")`
* `library(d3heatmap)`
* `install.packages("RColorBrewer")`
* `library(RColorBrewer)`
* `install.packages("htmlwidgets")`
* `library(htmlwidgets)`
* `install.packages("networkD3")`
* `library(networkD3)`
* `install.packages("data.table")`
* `library(data.table)`
* `install.packages("dplyr")`
* `library(dplyr)`
* `source("https://bioconductor.org/biocLite.R")`
* `biocLite("edgeR")`
* `library(edgeR)`
* `biocLite("GO.db")`
* `library(GO.db)`
* `biocLite("goseq")`
* `library(goseq)`
* `runApp("microscope")`

## Screenshots

<img width="1440" alt="heat" src="https://cloud.githubusercontent.com/assets/9893806/13304097/e36c1558-db20-11e5-8b79-83d69fd56ef0.png">

<img width="1440" alt="stats" src="https://cloud.githubusercontent.com/assets/9893806/13304104/ea4a1320-db20-11e5-9dbc-f4028f16c6c7.png">

<img width="1440" alt="go" src="https://cloud.githubusercontent.com/assets/9893806/13304107/eccd508a-db20-11e5-99c2-3ed8f24dd44a.png">

<img width="1440" alt="network" src="https://cloud.githubusercontent.com/assets/9893806/13304111/f0eb4d70-db20-11e5-83b1-b57c3f9a4985.png">
