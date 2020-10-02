# Code needed to install packages that are used to build the figures and results

install.packages("devtools")
devtools::install_github("Saeyslab/PeacoQC")
devtools::install_github("Saeyslab/FlowSOM", ref = "FlowSOM_v2")
devtools::install_github("jmeskas/flowCut")

install.packages("reshape2")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("tidyverse")
install.packages("cowplot")
install.packages("ggbeeswarm")
install.packages("viridis")
install.packages("ggpointdensity")
install.packages("xml2")
install.packages("ggrastr")
install.packages("svglite")
install.packages("pdist")


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("flowAI")
BiocManager::install("flowCore")
BiocManager::install("flowClean")
