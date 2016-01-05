# MicroScope

## About

MicroScope is an RShiny and JavaScript (D3.js) software program designed to produce dynamic, interactive heatmaps in a web browser. MicroScope allows you to magnify any portion of a heatmap by a simple click-and-drag feature to zoom in, and a click-once feature to zoom out. MicroScope is designed with large heatmaps in mind (e.g., gene expression heatmaps with thousands of genes), where individual entries quickly become unreadable as more and more add up. However, MicroScope allows you to repeatedly zoom in to any sector of the heatmap to investigate a region, cluster, or even a single gene. You can scroll up and down the page of your web browser to see more genes than fit your window. MicroScope also allows you to hover the mouse pointer over any specific gene to show precise expression level details.  In addition to visual magnification, MicroScope also allows users to analytically perform real-time statistical analyses with simple button clicks.  As such, MicroScope presents a significant advance in heatmap visualization technology over current standard protocols. 

MicroScope is an ongoing bioinformatics software project fully financially supported by the United States Department of Defense (DoD) through the National Defense Science and Engineering Graduate Fellowship (NDSEG) Program. This research was conducted with Government support under and awarded by DoD, Army Research Office (ARO), National Defense Science and Engineering Graduate (NDSEG) Fellowship, 32 CFR 168a.

Please cite: "Khomtchouk BB, Dargis-Robinson V, Hennessy JR, Wahlestedt C.  “MicroScope: real-time statistical analytics and visualization platform for gene expression heatmaps”. bioRxiv doi: http://dx.doi.org/10.1101/034694" within any source that makes use of any methods inspired by MicroScope. 

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
* `source("https://bioconductor.org/biocLite.R")`
* `biocLite("edgeR")`
* `library(edgeR)`
* `runApp("microscope")`

## Screenshots

<img width="1440" alt="picc1" src="https://cloud.githubusercontent.com/assets/9893806/12105623/74e1c930-b324-11e5-965e-08d523617f49.png">

<img width="1440" alt="picc2" src="https://cloud.githubusercontent.com/assets/9893806/12105624/771091be-b324-11e5-8676-e121ab2d84f7.png">

<img width="1440" alt="picc3" src="https://cloud.githubusercontent.com/assets/9893806/12105625/78933474-b324-11e5-8294-c6306f243f3b.png">
