library(PeacoQC)
library(flowCore)
library(flowAI)
library(flowClean)
library(flowCut)
library(tidyverse)
library(reshape2)

ff <- read.FCS("FlowCAP4/Compensated_transformed_fcs/010_comp_trans.fcs")

channels <- c(1,3,5:14, 18, 20, 21)

nrow(ff)

nr_cells_list <- c(1000,2500,5000,7500, 10000, 25000, 50000, 75000, 100000,
                   250000, 500000, 750000, 1000000, 1500000, 2000000, 2500000,
                   3000000)

ff@description$`$FIL` <- "010_comp_trans.fcs"


dir.create("Time_exp")

orig_expr_dataframe <- ff@exprs

while(nrow(ff) < 2500000){
  new_dataframe <- orig_expr_dataframe
  new_dataframe[,"Time"] <- orig_expr_dataframe[,"Time"] + 
    max(ff@exprs[,"Time"])
  ff@exprs <- rbind(ff@exprs, new_dataframe)
}


results <- list()

nr_cells <- 50000


for (run in c(1,2,3)){
  
  for (nr_cells in nr_cells_list ){
    
    print(nr_cells)
    
    subset_cells <- seq(1, nrow(ff), length.out = nr_cells)
    
    ff_subset <- ff[subset_cells,]
    
    time_peac <- system.time({peacoqc_res <- 
      PeacoQC(ff_subset, channels, output_directory = "Time_exp/",
              name_directory = paste0("PeacoQC_results", nr_cells),
              plot = F, save_fcs = F, report = F)})
    
    
    time_peac <- time_peac[3]
    
    out <- tryCatch({time_flowai <- 
      system.time(flowAI::flow_auto_qc(ff_subset,
                                       folder_results = paste0("Time_exp/FlowAI_",
                                                               nr_cells),
                                       ChExcludeFS = setdiff(colnames(ff),
                                                             colnames(ff)[channels]),
                                       html_report = F,
                                       fcs_QC = FALSE,
                                       mini_report = FALSE))},
      error = function(e){return(nr_cells)})
    
    
    if (class(out) == "numeric"){
      
      if (out == nr_cells){
        time_flowai <- NA} }else {
          time_flowai <- out[3]}
    
    out2 <- tryCatch({time_flowai_fs <- 
      system.time(flowAI::flow_auto_qc(ff_subset,
                                       folder_results = paste0("Time_exp/FlowAI_FS_", 
                                                               nr_cells),
                                       remove_from = "FS_FM",
                                       ChExcludeFS =  setdiff(colnames(ff),
                                                              colnames(ff)[channels]),
                                       html_report = F,
                                       fcs_QC = FALSE,
                                       mini_report = FALSE))},
      error = function(e){return(nr_cells)})
    
    
    if (class(out2) == "numeric"){
      if (out2 == nr_cells){
        time_flowai_FS <- NA} } else {
          time_flowai_FS <- out2[3]}
    
    
    time_flowcut <- system.time(flowCut(f = ff_subset, Channels = channels,
                                        Directory = paste0("Time_exp/flowCut_",nr_cells),
                                        Plot = "None"))
    
    time_flowcut <- time_flowcut[3]
    
    
    out3 <- tryCatch({time_flowclean <- 
      system.time(clean(ff_subset,
                        filePrefixWithDir = paste0("Time_exp/flowClean/",
                                                   ff@description$`$FIL`),
                        ext = ".fcs", vectMarkers = channels))},
      error = function(e){return(nr_cells)})
    
    if (class(out3) == "numeric"){
      if (class(out3) == nr_cells){
        time_flowclean <- NA}
    } else{
      time_flowclean <- time_flowclean[3]
    }
    
    results[[as.character(nr_cells)]] <- data.frame("Algorithm" = c("PeacoQC", 
                                                                    "flowAI", 
                                                                    "flowAI (FS & FM)", 
                                                                    "flowCut", 
                                                                    "flowClean"),
                                                    "Time" = c(time_peac, 
                                                               time_flowai,
                                                               time_flowai_FS, 
                                                               time_flowcut, 
                                                               time_flowclean),
                                                    "Nr_cells" = rep(nr_cells, 5))
    
  }
  
  results_matrix <- do.call(rbind, results)
  saveRDS(results_matrix, file = paste0("Time_exp/results_run",run,".RDS"))
  
}


result_files <- list.files("Time_exp/", pattern = "results", full.names = TRUE)


full_matrix <- do.call('cbind', lapply(result_files, readRDS))

results_matrix_new <- cbind(full_matrix[,c(1,3)], "Time" = apply(full_matrix[,c(2,5,8)], 1, mean))


results_matrix_new[which(results_matrix_new$Time <= 0.05), "Time"] <- NA
results_matrix_new[which(results_matrix_new$Time > 1000), "Time"] <- NA


p1 <- ggplot(data = results_matrix_new) +
  geom_point(aes(x = as.numeric(Nr_cells), y = as.numeric(Time, units = "secs"), 
                 colour = Algorithm, shape = Algorithm)) +
  geom_line(aes(x = as.numeric(Nr_cells), y = as.numeric(Time, units = "secs"), 
                colour = Algorithm)) +
  
  scale_x_log10(limits=c(1000, 3000000),
                breaks=c(1000, 10000, 100000, 1000000, 3000000),
                labels=c(1000, 10000, 100000, 1000000, 3000000)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks=c(0,5,10,50,250,750),
                     labels=c(0,5,10,50,250,750)) + 
  scale_shape_manual(values=c(19,2,15,3,11))+
  theme_minimal() +
  labs(x = "Nr of events", y = "Time (s)")
p1



ggsave(plot = p1, 
       filename = "Time_exp/Figure7.png",
       width = 8, height = 3)













