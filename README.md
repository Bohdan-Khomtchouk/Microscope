# MicroScope

## About

MicroScope is an RShiny and JavaScript (D3.js) software program designed to produce dynamic, interactive heatmaps in a web browser. MicroScope allows you to magnify any portion of a heatmap by a simple click-and-drag feature to zoom in, and a click-once feature to zoom out. MicroScope is designed with large heatmaps in mind (e.g., gene expression heatmaps with thousands of genes), where individual entries quickly become unreadable as more are added. However, MicroScope allows you to repeatedly zoom in to any sector of the heatmap to investigate a region, cluster, or even a single gene. You can scroll up and down the page of your web browser to see more genes than fit your window. MicroScope also allows you to hover the mouse pointer over any specific gene to show precise expression level details. As such, MicroScope presents a significant advance in heatmap visualization technology over current standard protocols.       

## Installation

### Requirements

* R programming language
  * RStudio

## How to run

##### Git clone this repo to your computer, and in RStudio type:
* `setwd("~/path/to/my_directory/MICROSCOPE")`
* `install.packages("d3heatmap")`
* `install.packages("shiny")`
* `library(shiny)`
* Now go to ui.r and manually change `pwd<-"/path/to/my_directory/MICROSCOPE/genes_file.csv"` to your own file and path
* `runApp("microscope")`

## Screenshots

<img width="1440" alt="image1" src="https://cloud.githubusercontent.com/assets/9893806/11605850/75c634f8-9ad8-11e5-83af-3cfa6a6387e5.png">

<img width="1440" alt="image2" src="https://cloud.githubusercontent.com/assets/9893806/11605851/7b3ddae4-9ad8-11e5-8f15-51e16b581beb.png">

<img width="1440" alt="image3" src="https://cloud.githubusercontent.com/assets/9893806/11605852/7ce2319c-9ad8-11e5-8d40-e2dfadd767f0.png">

<img width="1264" alt="image4" src="https://cloud.githubusercontent.com/assets/9893806/11605854/7ed04dea-9ad8-11e5-8556-3f76e9ee63b4.png">

<img width="1440" alt="image5" src="https://cloud.githubusercontent.com/assets/9893806/11605855/8088db8e-9ad8-11e5-90fe-651806b87eb4.png">

<img width="1439" alt="image6" src="https://cloud.githubusercontent.com/assets/9893806/11605857/8217f5e8-9ad8-11e5-8187-4c052254b684.png">
