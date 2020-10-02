library(PeacoQC)
library(flowClean)
library(flowCut)
library(flowAI)
library(flowCore)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(reshape2)

# Run different QC algorithms

files <- list.files("FlowCAP4/Compensated_transformed_fcs/")


channels <- c(1,3,5,7,11,14)

file <- files[2]

rename <- read.table("FlowCAP4/SVG_files_manual_annotations/FlowCAP_svgNames.txt", 
                     header = F, 
                     sep = "")[,c(1,3)]
rename[,1] <- sprintf("%03d", rename[,1]) 

files <- files[sub("_comp_trans","",files)%in%  paste0(rename[,1], ".fcs")]
file <- files[54]



for (file in files[1:55]){
    ff <- read.FCS(paste0("FlowCAP4/Compensated_transformed_fcs/", file))
    
    # PeacoQC analyse
    
    PeacoQC_res <- PeacoQC(ff,
                           channels,
                           output_directory = "FlowCAP4",
                           name_directory = "PeacoQC_Results",
                           plot = T)

    # flowclean analyse

    # out <- tryCatch({flowclean_res <- flowClean::clean(ff, channels,
    #                                   filePrefixWithDir = paste0("FlowCAP4/flowClean/",
    #                                                              sub(".fcs", "",file)),
    #                                   ext = "fcs", diagnostic = TRUE)
    # 
    # write.FCS(flowclean_res, paste0("FlowCAP4/flowClean/", file))}
    # )
    # 
    # # flowCut analyse
    # 
    # flowcut_res <- flowCut(ff, Channels = channels, 
    #                        Directory = "FlowCAP4/flowCut",Plot = "All")
    # 
    # saveRDS(flowcut_res, file = paste0("FlowCAP4/flowCut/",
    #                                        sub(".fcs", ".Rds", file)))
    # 
    # write.FCS(flowcut_res$frame,
    #           filename = paste0("FlowCAP4/flowCut/", basename(file)))
    # 
    # 
    # # FlowAI one step
    # 
    # # Make sure that filename is stored correctly
    # 
    # ff@description$`$FIL` <- file
    # 
    # flowai_full <- flow_auto_qc(ff,
    #                             ChExcludeFM = colnames(ff)[-channels],
    #                             ChExcludeFS = colnames(ff)[-channels],
    #                             folder_results = "FlowCAP4/FlowAI_full",
    #                             output = 2)
    

    
}
