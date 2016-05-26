# MicroScope

## About

MicroScope is a user-friendly ChIP-seq and RNA-seq software suite designed for the interactive visualization and analysis of gene expression heatmaps, including integrated features to support: principal component analysis, differential expression analysis, gene ontology analysis, and dynamic network visualization of genes directly from a heatmap.

MicroScope is financially supported by the United States Department of Defense (DoD) through the National Defense Science and Engineering Graduate Fellowship (NDSEG) Program. This research was conducted with Government support under and awarded by DoD, Army Research Office (ARO), National Defense Science and Engineering Graduate (NDSEG) Fellowship, 32 CFR 168a.

Please cite: "Khomtchouk BB, Hennessy JR, Dargis-Robinson V, Wahlestedt C.  “MicroScope: ChIP-seq and RNA-seq software analysis suite for gene expression heatmaps” (submitted). bioRxiv doi: http://dx.doi.org/10.1101/034694" within any source that makes use of any methods inspired by MicroScope. 

## Usage (for general public)

##### Just click here!... http://microscopebioinformatics.org/

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
* `biocLite("org.Bt.eg.db")`
* `biocLite("org.Ce.eg.db")`
* `biocLite("org.Cf.eg.db")`
* `biocLite("org.Dm.eg.db")`
* `biocLite("org.Dr.eg.db")`
* `biocLite("org.Hs.eg.db")`
* `biocLite("org.Mm.eg.db")`
* `biocLite("org.Pt.eg.db")`
* `biocLite("org.Rn.eg.db")`
* `biocLite("org.Ss.eg.db")`
* `biocLite("org.Sc.sgd.db")`
* `runApp("microscope")`

## Screenshots

<img width="1440" alt="heat" src="https://cloud.githubusercontent.com/assets/9893806/13304097/e36c1558-db20-11e5-8b79-83d69fd56ef0.png">

<img width="1440" alt="stats" src="https://cloud.githubusercontent.com/assets/9893806/13304104/ea4a1320-db20-11e5-9dbc-f4028f16c6c7.png">

<img width="1440" alt="go" src="https://cloud.githubusercontent.com/assets/9893806/13304107/eccd508a-db20-11e5-99c2-3ed8f24dd44a.png">

<img width="1440" alt="network" src="https://cloud.githubusercontent.com/assets/9893806/13304111/f0eb4d70-db20-11e5-83b1-b57c3f9a4985.png">
