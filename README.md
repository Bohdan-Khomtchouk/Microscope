![microscope_logo](https://cloud.githubusercontent.com/assets/9893806/16903838/2b0f7b82-4c57-11e6-875a-dbe8f4f5ae11.png)

# MicroScope

## About

We propose a user-friendly ChIP-seq and RNA-seq software suite for the interactive visualization and analysis of genomic data, including integrated features to support differential expression analysis, interactive heatmap production, principal component analysis, gene ontology analysis, and dynamic network visualization.

MicroScope is financially supported by the United States Department of Defense (DoD) through the National Defense Science and Engineering Graduate Fellowship (NDSEG) Program. This research was conducted with Government support under and awarded by DoD, Army Research Office (ARO), National Defense Science and Engineering Graduate (NDSEG) Fellowship, 32 CFR 168a.

Please cite: "Khomtchouk BB, Hennessy JR, Wahlestedt C.  “MicroScope: ChIP-seq and RNA-seq software analysis suite for gene expression heatmaps” BMC Bioinformatics.  2016 (in press). bioRxiv doi: http://dx.doi.org/10.1101/034694" within any source that makes use of any methods inspired by MicroScope. 

## Usage (for general public)

##### http://microscopebioinformatics.org/

## Installation (for developers only)

### Requirements for developers

* R programming language
  * RStudio

## How to run (for developers only)

##### Git clone this repo to your computer, and in RStudio type:
* `setwd("~/path/to/my_directory/Microscope")`
* `install.packages("shiny")`
* `library(shiny)`
* `library(shinyapps)`
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

![fig1](https://cloud.githubusercontent.com/assets/9893806/16438422/aac3ba7a-3d7d-11e6-867d-8aaf8d6b7e86.png)
<img width="655" alt="fig1 1" src="https://cloud.githubusercontent.com/assets/9893806/16439227/4d076060-3d84-11e6-9d39-4d54728ed232.png">
![fig2](https://cloud.githubusercontent.com/assets/9893806/16438431/bc6c3a68-3d7d-11e6-9fdf-e570865c906d.png)
![fig3](https://cloud.githubusercontent.com/assets/9893806/16438436/c1a096fa-3d7d-11e6-9c35-d80c97dd0e2d.png)
![fig4](https://cloud.githubusercontent.com/assets/9893806/16438438/c4be0be2-3d7d-11e6-9381-1101206ade6d.png)
![fig5](https://cloud.githubusercontent.com/assets/9893806/16438441/c999fd10-3d7d-11e6-8dd6-65d37ab43e06.png)
![fig6](https://cloud.githubusercontent.com/assets/9893806/16438443/ceb93b62-3d7d-11e6-817c-c12dee745c9b.png)
