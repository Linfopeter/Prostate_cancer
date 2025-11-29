#  Paquetes necesarios para análisis de RNAseq
#  ===========================================================

# Install Rtools
# https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html

install.packages("pkgbuild")
pkgbuild::check_build_tools(debug = TRUE)

install.packages("tidyverse")
install.packages("dplyr")
install.packages("plotly")
install.packages("enrichR")
install.packages("ashr")
install.packages("VennDiagram")
install.packages("gplots")
install.packages("pheatmap")
install.packages("WGCNA")
install.packages("gridExtra")
install.packages("igraph")
install.packages("devtools")

# Instalar paquetes de Bioconductor
# ===========================================================

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("BiocVersion")
BiocManager::install("ReportingTools")
#BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
#BiocManager::install("Rsubread")
#BiocManager::install("S4Vectors")
#BiocManager::install("DESeq2")
#BiocManager::install("gage")
#BiocManager::install("biomaRt")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("GO.db")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install('EnhancedVolcano')
#BiocManager::install("apeglm")
#BiocManager::install("vsn")
#BiocManager::install("enrichplot")
#BiocManager::install("impute")
#BiocManager::install("GEOquery")
#BiocManager::install('PCAtools')
#BiocManager::install("GWENA")
#BiocManager::install('mixOmics')
#BiocManager::install('igraph')
#BiocManager::install("fields")
#BiocManager::install("dynamicTreeCut")
#BiocManager::install("qvalue")
#BiocManager::install("flashClust")
#BiocManager::install("Hmisc")
#BiocManager::install("WGCNA")
#BiocManager::install("pathview")
#BiocManager::install("RCy3")
#BiocManager::install("Rfastp")



#Para actualizar los paquetes de Bioconductor
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install()

# Instalacion de paquetes Devtools
# =============================================================
#if (!require(devtools))
#install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
devtools::install_github("https://github.com/hancockinformatics/pathlinkR")
devtools::install_github("noriakis/ggkegg")

# Instalar PathfindR (analisis de enriquecimiento)
# =============================================================
install.packages("pak") # if you have not installed "pak"
pak::pkg_install("pathfindR")

remotes::install_github("kevinblighe/CorLevelPlot")



