# MicroScope

## About

MicroScope is a user-friendly ChIP-seq and RNA-seq software suite designed for the interactive visualization and analysis of gene expression heatmaps, including integrated features to support: principal component analysis, differential expression analysis, gene ontology analysis, and dynamic network visualization of genes directly from a heatmap.

MicroScope is financially supported by the United States Department of Defense (DoD) through the National Defense Science and Engineering Graduate Fellowship (NDSEG) Program. This research was conducted with Government support under and awarded by DoD, Army Research Office (ARO), National Defense Science and Engineering Graduate (NDSEG) Fellowship, 32 CFR 168a.

Please cite: "Khomtchouk BB, Hennessy JR, Dargis-Robinson V, Wahlestedt C.  “MicroScope: ChIP-seq and RNA-seq software analysis suite for gene expression heatmaps” (submitted). bioRxiv doi: http://dx.doi.org/10.1101/034694" within any source that makes use of any methods inspired by MicroScope. 

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

![fig1](https://cloud.githubusercontent.com/assets/9893806/15684804/cde007c8-2736-11e6-9893-5925d5e15561.png)

![fig2](https://cloud.githubusercontent.com/assets/9893806/15684822/f82654c4-2736-11e6-881c-7c4aa47a1960.png)

![fig3](https://cloud.githubusercontent.com/assets/9893806/15684826/fc0498c6-2736-11e6-94a5-358146b825a8.png)

![fig4](https://cloud.githubusercontent.com/assets/9893806/15684834/04ff3f76-2737-11e6-8c67-d68da29fd676.png)

![fig5](https://cloud.githubusercontent.com/assets/9893806/15684844/0d573c78-2737-11e6-9c0e-c432a8408b3d.png)

![fig6](https://cloud.githubusercontent.com/assets/9893806/15684852/1620700e-2737-11e6-8f06-2c162e247983.png)
