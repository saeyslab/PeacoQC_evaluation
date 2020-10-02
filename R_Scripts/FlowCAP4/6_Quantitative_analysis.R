library(flowCore)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggbeeswarm)
library(ggpubr)

Se_SP <- function(kept, trues, falses){
    sensitivity <- sum(falses[!kept])/(
        sum(falses[!kept]) + sum(falses[kept])) 
    
    specificity <- sum(trues[kept])/( 
        sum(trues[!kept]) + sum(trues[kept]))

    
    if(is.nan(sensitivity)){
        sensitivity <- 1
    }
    
    return(list("sensitivity" = sensitivity, "specificity" = specificity))
    
}



data_total <- list()
files <- list.files("FlowCAP4/Compensated_transformed_fcs/")

channels <- c(1,3,5,7,11,14)


rename <- read.table("FlowCAP4/SVG_files_manual_annotations/FlowCAP_svgNames.txt", 
                     header = F, 
                     sep = "")[,c(1,3)]
rename[,1] <- sprintf("%03d", rename[,1]) 

files <- files[sub("_comp_trans","",files)%in%  paste0(rename[,1], ".fcs")]

file <- files[44]


for (file in files){
    
    print(file)
    
    ff <- read.FCS(file.path("FlowCAP4/Compensated_transformed_fcs/", file))
    
    set.seed(1)
    
    sampled_cells <- sample(nrow(ff), 10000)

    table_file <- rename[grep(sub("_comp_trans.fcs","",file), rename[,1]),2]
    
    table_file <- readRDS(paste0("FlowCAP4/Full_annotated_tables/", 
                                 table_file, ".Rds"))
    
    
    total_table <- matrix(nrow = nrow(ff), ncol = ncol(table_file), NA)
    total_table[,1] <- 1:nrow(total_table)
    
     annotator <- 2
    for(annotator in 2:ncol(table_file)){
        total_table[sort(sampled_cells),annotator] <- pull(table_file, annotator)
        total_table[,annotator] <- zoo::na.locf(total_table[,annotator], na.rm = F)
        total_table[,annotator] <- zoo::na.locf(total_table[,annotator], 
                                           fromLast = T, na.rm = F)
    }
    
    
    total_values <- rowSums(total_table[,-1])
    
    false_values <- rep(0,nrow(total_table))
    
    false_values[which(total_values > (ncol(total_table) - 1)/2)] <- 1
    
    true_values <- 1 - false_values
    
    random_cells <- rep(TRUE, nrow(ff))
    
    random_cells[sample(nrow(ff), length(which(false_values ==1)))] <- FALSE
    
    flowai_ff <- read.FCS(paste0("FlowCAP4/FlowAI_full/",
                                 sub(".fcs", "_QC.fcs", file)))
    flowai_kept <- flowai_ff@exprs[,"remove_from_all"] < 10000
    
    flowai_ff <- flowai_ff[which(flowai_ff@exprs[,"remove_from_all"] < 10000), ]
    
    peacoqc_ff <- read.FCS(paste0("FlowCAP4/PeacoQC_Results/fcs_files/", 
                                  sub(".fcs","_QC.fcs",file)))
    
    peacoqc_kept <- rep(FALSE, max(ff@exprs[,"Original_ID"]))
    peacoqc_kept[peacoqc_ff@exprs[,"Original_ID"]] <- TRUE
    peacoqc_kept <- peacoqc_kept[ff@exprs[,"Original_ID"]]
    
    
    flowaifs_ff <- read.FCS(paste0("FlowCAP4/FlowAI_FS_FM/",
                                   sub(".fcs","_QC.fcs",file)))
    flowaifs_kept <- flowaifs_ff@exprs[,"remove_from_FS_FM"] < 10000
    
    flowaifs_ff <- flowaifs_ff[which(flowaifs_ff@exprs[,"remove_from_FS_FM"] <
                                         10000), ]
    
    flowcut_kept <- rep(TRUE, nrow(ff))
    flowcut_res <- readRDS(paste0("FlowCAP4/flowCut/",
                                  sub(".fcs", ".Rds", file)))
    flowcut_kept[flowcut_res$ind] <- FALSE
    
    flowclean_ff <- read.FCS(paste0("FlowCAP4/flowClean/",file))
    flowclean_ketp <- flowclean_ff@exprs[,"GoodVsBad"] < 10000

    random_SP <- Se_SP(random_cells, true_values, false_values)
    peacoqc_SP <- Se_SP(peacoqc_kept, true_values, false_values)
    flowai_SP <- Se_SP(flowai_kept, true_values, false_values)
    flowaifs_SP <- Se_SP(flowaifs_kept, true_values, false_values)
    flowcut_SP <- Se_SP(flowcut_kept, true_values, false_values)
    flowclean_SP <- Se_SP(flowclean_ketp, true_values, false_values)
    
    
    cells_random <- ((nrow(ff) - length(which(random_cells == TRUE)))/nrow(ff)) *100
    cells_FlowAI <- ((nrow(ff) - nrow(flowai_ff))/nrow(ff))*100
    cells_peacoqc <-((nrow(ff) - nrow(peacoqc_ff))/nrow(ff))*100
    cells_FlowAIFS <- ((nrow(ff) - nrow(flowaifs_ff))/nrow(ff))*100
    cells_flowcut <-((nrow(ff) - nrow(flowcut_res$frame))/nrow(ff))*100
    cells_flowclean <- ((nrow(ff) -
                             nrow(flowclean_ff[flowclean_ff@exprs[,"GoodVsBad"] < 10000,]))/
                            nrow(ff)) *100
    
    cells_original <- apply(table_file[,-1], 2,
                            function(x)((nrow(ff) - nrow(ff[!x]))/nrow(ff))*100)
    
    
    overview_data <- data.frame(
        Algorithm = c("Random (Percentage removed = golden standard)",
                      "PeacoQC", "flowAI", "flowAI (FS & FM)", 
                      "flowCut", "flowClean",
                      rep("Original", ncol(table_file) - 1)),
        Sensitivity = c(random_SP$sensitivity, peacoqc_SP$sensitivity, 
                        flowai_SP$sensitivity,
                        flowaifs_SP$sensitivity, flowcut_SP$sensitivity, 
                        flowclean_SP$sensitivity,
                        rep(1, ncol(table_file) -1)),
        Specificity = c(random_SP$specificity, peacoqc_SP$specificity, 
                        flowai_SP$specificity,
                        flowaifs_SP$specificity, flowcut_SP$specificity, 
                        flowclean_SP$specificity,
                        rep(1, ncol(table_file) -1)),
        Filename = rep(file, 6 + ncol(table_file) -1),
        perc_removed = c(cells_random, cells_peacoqc, cells_FlowAI, 
                         cells_FlowAIFS, cells_flowcut,
                         cells_flowclean, cells_original))
    
    
    data_total[[file]] <- overview_data

    
}




# --------------------------- Overview plot -----------------------------------

overview_data <- do.call(rbind, data_total) %>% 
    replace_na(list(CV = 0)) %>%
    mutate(Filename = factor(Filename)) %>% 
    filter(Algorithm != "Original") %>% 
    mutate("Percentage Kept" = (1 - (perc_removed/100))) %>% 
    mutate("Balanced accuracy" = (Sensitivity + Specificity)/2) %>% 
    rename( `Percentage Removed` = perc_removed) %>% 
    mutate("Percentage Removed" = `Percentage Removed`/100)

overview_data_test <- melt(overview_data, 
                           measure.vars = c("Percentage Removed", 
                                            "Sensitivity", 
                                            "Specificity", 
                                            "Balanced accuracy"))


overview_data_test$variable <- factor(overview_data_test$variable, 
                                      levels = rev(levels(overview_data_test$variable)))



medians_figure <- overview_data %>% group_by(Algorithm) %>% 
    summarise(median(`Balanced accuracy`),
              median(Sensitivity),
              median(Specificity),
              median(`Percentage Removed`))
colnames(medians_figure) <- c("Algorithm",
                              "Balanced accuracy",
                              "Sensitivity",
                              "Specificity",
                              "Percentage Removed")


medians_figure <- medians_figure %>% 
    gather(key = "Measurements", value = "value", -"Algorithm")

medians_dataframe <- as.data.frame(medians_figure)

medians_dataframe$value <- medians_dataframe$value %>% round(.,3) %>% format(.,nsmall = 2)


p <- ggplot(overview_data_test, aes(x = variable, y = value))+
     geom_violin(scale = "width", aes(fill = variable)) + 
    geom_quasirandom(size = 0.2, color = "grey30") +
    theme_bw()+ 
    facet_wrap(~Algorithm, ncol = 1) + 
    theme(axis.title.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          strip.text = element_text(hjust = 0, face = "bold"),
          strip.background = element_rect(fill = "white", colour = "#00000000"),
          axis.text.y = element_text(margin = margin(r = 40))) +
    coord_flip(clip = "off", ylim = c(0,1)) + 
    geom_text(data = medians_dataframe, aes(y = -0.13, x = Measurements, 
                                            label = as.character(value)), 
              size = 3)



ggsave(plot =p, "Figure6.png", height = 10, width = 6)


# ------------------------------- Statistical analysis ------------------------


different_statistics <- c("A. Balanced accuracy", "B. Sensitivity", 
                          "C. Specificity", "D. Percentage Removed")

plot_list <- list()

statistic <- different_statistics[1]

for (statistic in different_statistics){
    
    df <- overview_data_test %>% filter(variable == sub(".\\. ","",statistic))
    
    
    pwt <- pairwise.wilcox.test(df$value, 
                                df$Algorithm,
                                p.adjust.method = "BH")
    
    m <- formatC(pwt$p.value, format = "e", digits = 2)
    m[which(m == " NA")] <- ""
    rownames(m)[5] <- "Random"
    
    
    p <- ggtexttable(m, 
                     theme = ttheme(base_size = 10,
                                    rownames.style = 
                                        colnames_style(face = "bold", size = 10),
                                    tbody.style = 
                                        tbody_style(fill = "white", size = 10)))

    
        
    for(row in 1:5){
        for(col_value in 1:5){
            val <- m[row, col_value]
            if(val != "" & as.numeric(val) < 0.05){
                p <- table_cell_bg(p, row = row+1, 
                                   column = col_value+1, 
                                   color = "white", 
                                   fill = "#FF000022")
            }
        }
    }
    
    plot_list[[statistic]] <- p  %>% 
        tab_add_title(text = statistic, face = "bold")
    
    
    
    
}

final_plot <- ggarrange(plotlist = plot_list, nrow = 2, ncol = 2,
          align = "hv")

ggsave(plot = final_plot, "Supp_figure_2.png", width = 12, height = 6)



