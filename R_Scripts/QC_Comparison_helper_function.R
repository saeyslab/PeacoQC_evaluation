# ---------------------------- Extract cells ----------------------------------

ExtractQCCells <- function(file, ff, directory){
    
    flowai_ff <- tryCatch({read.FCS(paste0(directory,"/FlowAI_full/",
                                           sub(".fcs", "_QC.fcs", file)))},
                          error = function(e){return(NA)})
    
    if (is(flowai_ff, "flowFrame")){
        flowai_kept <- flowai_ff@exprs[,"remove_from_all"] < 10000
        flowai_ff <- flowai_ff[which(flowai_ff@exprs[,"remove_from_all"] < 10000), ]
        cells_flowAI <- ((nrow(ff) - nrow(flowai_ff))/nrow(ff))*100
        
    } else{ flowai_kept <- rep(FALSE, nrow(ff))
    }
    
    
    peacoqc_ff <- read.FCS(paste0(directory,"/PeacoQC_Results/fcs_files/", 
                                  sub(".fcs","_QC.fcs",file)))
    
    peacoqc_kept <- rep(FALSE, max(ff@exprs[,"Original_ID"]))
    peacoqc_kept[peacoqc_ff@exprs[,"Original_ID"]] <- TRUE
    peacoqc_kept <- peacoqc_kept[ff@exprs[,"Original_ID"]]
    
    flowcut_kept <- rep(FALSE, nrow(ff))
    flowcut_ff <- read.FCS(paste0(directory,"/flowCut/", sub(".fcs", "_QC.fcs", file)))
    flowcut_kept[which(ff@exprs[,"Original_ID"] %in% flowcut_ff@exprs[,"Original_ID"])] <- TRUE
    
    flowclean_ff <- read.FCS(paste0(directory,"/flowClean/",sub(".fcs", "_QC.fcs", file)))
    flowclean_ketp <- flowclean_ff@exprs[,"GoodVsBad"] < 10000
    
    frame_cells <- data.frame("idc" = seq_len(nrow(ff)),
                              "PeacoQC" = peacoqc_kept,
                              "flowCut" = flowcut_kept,
                              "flowAI" = flowai_kept,
                              "flowClean" = flowclean_ketp)
    
    return(frame_cells)
}


# ------------------------------- Calculate median and quantiles --------------

CalculateMedianQuantile <- function(ff, channel){
    
    if (ff@description$FILENAME == "Data_flowrepository/FR-FCM-Z2EY/PP/Full staining A1_PP.fcs"){
        
        breaks <- cut(seq_along(ff@exprs[,channel]),
                      breaks = seq(0, nrow(ff@exprs)+100, by = 100),
                      labels = FALSE)
        mid_breaks <- seq(50, 
                          length(ff@exprs[,channel]) + 50, by = 100)
    } else{
        
        breaks <- cut(seq_along(ff@exprs[,channel]),
                      breaks = seq(0, nrow(ff@exprs)+1000, by = 1000),
                      labels = FALSE)
        mid_breaks <- seq(500, 
                          length(ff@exprs[,channel]) + 500, by = 1000)
    }
    
    medians <- tapply(ff@exprs[,channel], breaks, median)
    q25 <- tapply(ff@exprs[,channel], breaks, quantile, 0.25)
    q75 <- tapply(ff@exprs[,channel], breaks, quantile, 0.75)
    
    mid_breaks <- mid_breaks[1:length(medians)]
    
    
    median_quantile_frame <- data.frame("idc" = mid_breaks,
                                        "median" = medians,
                                        "q25" = q25,
                                        "q75" = q75)
    
    return(median_quantile_frame)
    
}


# ----------------- Heatmap QC -------------------------------------------------

HeatmapQC <- function(frame_cells.melted){
    p2 <- ggplot()+
        geom_tile(data = frame_cells.melted,aes(x = idc, 
                                                y = variable, 
                                                fill = value), 
                  height = 0.9) +
        theme(legend.position="bottom") +
        theme_minimal() +
        xlab("Events ordered on time") +
        theme(axis.title.y = element_blank(),
              legend.title = element_blank(),
              axis.text.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(size = 13)) +
        scale_fill_manual(labels = c("FALSE" = "Removed by algorithm",
                                     "TRUE" = "Kept by algorithm",
                                     "NA" = "Algorithm did not work"),
                          values = c("TRUE" = "#27476E",
                                     "FALSE" = "#E26D5C",
                                     "NA" = "#DF2935")) 
    return(p2)
}


# --------------------------- Heatmap mass cytometry ---------------------------
HeatmapQCMassData <- function(frame_cells.melted){
    frame_cells.melted[is.na(frame_cells.melted)] <- "NA"
    
    frame_cells.melted$value <- factor(frame_cells.melted$value, levels = 
                                           c( "TRUE", "FALSE", "NA"))
    
    p2 <- ggplot()+
        geom_tile(data = frame_cells.melted,aes(x = idc, 
                                                y = variable, 
                                                fill = value), 
                  height = 0.9) +
        theme(legend.position="bottom") +
        theme_minimal() +
        xlab("Events ordered on time") +
        theme(axis.title.y = element_blank(),
              legend.title = element_blank(),
              axis.text.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(size = 13)) +
        scale_fill_manual(labels = c("FALSE" = "Removed by algorithm",
                                     "TRUE" = "Kept by algorithm",
                                     "NA" = "Algorithm did not work"),
                          values = c("TRUE" = "#27476E",
                                     "FALSE" = "#E26D5C",
                                     "NA" = "grey55")) 
    return(p2)
}

# -------------------------- Build background blocks ---------------------------

BuildBackgroundBlocks <- function(h, test, file){
    
    
    if (file == "flow_file_afterdapi.fcs"){
        run_length <- rle(rep("Normal" , length(h$counts)))
        fill_blocks <- run_length$values
        
    } else if (file == "flow_file_low_medium_high_speed.fcs"){
        
        perturbations <- ifelse(h$counts > 2000 , "Perturbation 2",
                                ifelse(h$counts > 1000, "Perturbation 1",
                                       "Normal"))
        run_length <- rle(perturbations)
        fill_blocks <- run_length$values
        
        
    } else if (file ==  "flow_file_high_low_high_speed.fcs"){
        run_length <- rle(h$counts > 500)
        fill_blocks <- ifelse(run_length$values == TRUE, "Normal",
                              "Perturbation")
        
    } else if( file ==  "flow_file_start_up.fcs"){
        run_length <- rle(c(rep(FALSE, 2), rep(TRUE, 14)))
        fill_blocks <- ifelse(run_length$values == TRUE, "Normal",
                              "Perturbation")
        
    } else if (file == "Mass_file.fcs"){
        run_length <- rle(rep("Normal" , length(h$counts)))
        fill_blocks <- run_length$values
    } else if (file == "Spectral_file.fcs"){
        run_length <- rle(rep("Normal" , length(h$counts)))
        fill_blocks <- run_length$values
    }
    
    
    x_min <- test[c(1, cumsum(run_length$lengths)
                    [-length(run_length$lengths)] + 1)]
    x_max <- test[cumsum(run_length$lengths) +1]
    
    
    overview_blocks <- data.frame(x_min=x_min,
                                  x_max=x_max,
                                  y_min=-Inf,
                                  y_max=Inf,
                                  fill_blocks=fill_blocks)
    
    blocks_plot <- geom_rect(data=overview_blocks,
                             mapping=aes(xmin=x_min,
                                         xmax=x_max,
                                         ymin=y_min,
                                         ymax=y_max,
                                         fill=fill_blocks),
                             alpha=0.4,
                             show.legend=TRUE)
    
    return(blocks_plot)
}

# --------------------------- BuildManualScale --------------------------------

BuildManualScale <- function(file){
    
    if (file %in% c("flow_file_afterdapi.fcs",
                    "Mass_file.fcs",
                    "Spectral_file.fcs")){
        scale_plot <- scale_fill_manual(name="",
                                        limits = "",
                                        values=c(`Normal`="white"),
                                        guide=guide_legend(override.aes=
                                                               list(alpha=0.4)))
        
    } else if (file == "flow_file_low_medium_high_speed.fcs"){
        scale_plot <- scale_fill_manual(name="",
                                        limits = c("Perturbation 1", "Perturbation 2"),
                                        values=c(`Perturbation 1`="#A9B4C2",
                                                 `Perturbation 2` = "grey48", 
                                                 `Normal`="white"),
                                        guide=guide_legend(override.aes=
                                                               list(alpha=0.4)))
        
        
    } else {
        scale_plot <- scale_fill_manual(name="",
                                        limits = "Perturbation",
                                        values=c(`Perturbation`="#A9B4C2",
                                                 `Normal`="white"),
                                        guide=guide_legend(override.aes=
                                                               list(alpha=0.4)))
        
        
    }
    
    return(scale_plot)
    
    
}



MakeTitle <- function(file){
    if (file == "flow_file_afterdapi.fcs"){
        title <- "Measured after DAPI staining without cleaning"
        
    } else if (file == "flow_file_low_medium_high_speed.fcs"){
        
        title <- "Speed change (low-median-high)"
        
    } else if (file ==  "flow_file_high_low_high_speed.fcs"){
        title <- "Speed change (high-low-high-low)"
        
    } else if( file ==  "flow_file_start_up.fcs"){
        title <- "Start up of sample"
        
    } else if (file == "Mass_file.fcs"){
        title <- "Mass cytometry data"
        
    } else if (file == "Spectral_file.fcs"){
        title <- "Spectral cytometry data"
    }
    
    return(title)
    
}


RemoveMargins_masscytometry <- function(
    ff,
    channels,
    channel_specifications=NULL,
    output="frame") {
    
    
    meta <- flowWorkspace::pData(flowCore::parameters(ff))
    rownames(meta) <- meta[, "name"]
    
    if(!is.null(channel_specifications)){
        meta[names(channel_specifications),
             c("minRange", "maxRange")] <- do.call(rbind, channel_specifications)
    }
    
    # Initialize variables
    selection <- rep(TRUE, times=dim(ff)[1])
    e <- flowCore::exprs(ff)
    
    if(is.numeric(channels)){
        channels <- colnames(flowCore::exprs(ff))[channels]
    }
    # Make selection
    for (d in channels) {
        
        selection <- selection &
            e[, d] >= max(min(meta[d, "minRange"], 0), min(e[, d])) &
            e[, d] < min(min(meta[d, "maxRange"], 6.5), max(e[, d]))
    }
    
    
    if (length(which(selection == FALSE))/length(selection) > 0.1) {
        warning(PeacoQC:::StrMessage(c("More then ",
                                       round(length(which(selection == FALSE))/length(selection) * 100, 2),
                                       "% is considered as a margin event in file ",
                                       basename(flowCore::keyword(ff)$FILENAME),
                                       ". This should be verified.")))
    }
    
    new_ff <- ff[selection, ]
    new_ff <- PeacoQC:::AppendCellID(new_ff, which(selection))
    
    if (output == "full"){
        return(
            list("flowframe"=new_ff,
                 "indices_margins"=which(selection == FALSE)))
    } else if (output == "frame"){
        return(new_ff)
    }
}


ExtractQCCellsMass <- function(file, ff, directory){
    
    flowai_ff <- tryCatch({read.FCS(paste0(directory, "/FlowAI_res/",
                                           sub(".fcs", "_QC.fcs", file)))},
                          error = function(e){return(NA)})
    
    if (is(flowai_ff, "flowFrame")){
        flowai_kept <- flowai_ff@exprs[,"remove_from_all"] < 10000
        flowai_ff <- flowai_ff[which(flowai_ff@exprs[,"remove_from_all"] < 10000), ]
        cells_flowAI <- ((nrow(ff) - nrow(flowai_ff))/nrow(ff))*100
        
    } else{ flowai_kept <- rep(NA, nrow(ff))
    }
    
    
    peacoqc_ff <- read.FCS(paste0(directory, "/PeacoQC_Results/fcs_files/", 
                                  sub(".fcs","_QC.fcs",file)))
    
    peacoqc_kept <- rep(FALSE, max(ff@exprs[,"Original_ID"]))
    peacoqc_kept[peacoqc_ff@exprs[,"Original_ID"]] <- TRUE
    peacoqc_kept <- peacoqc_kept[ff@exprs[,"Original_ID"]]
    
    flowcut_ff <- tryCatch({read.FCS(paste0(directory, "/flowCut/", file))},
                           error = function(e){return(NA)})
    if (!is.na(flowcut_ff)){
        flowcut_kept <- rep(FALSE, nrow(ff))
        
        flowcut_kept[which(ff@exprs[,"Original_ID"] %in% flowcut_ff@exprs[,"Original_ID"])] <- TRUE
    } else{flowcut_kept <- rep(NA, nrow(ff))}
    
    
    flowclean_ff <- tryCatch({read.FCS(paste0(directory, "/flowClean/",file))},
                             error = function(e){return(NA)})
    
    if (is(flowclean_ff, "flowFrame")){
        flowclean_ketp <- flowclean_ff@exprs[,"GoodVsBad"] < 10000
        
        
    } else{ flowclean_ketp <- rep(NA, nrow(ff))
    }
    
    
    frame_cells <- data.frame("idc" = seq_len(nrow(ff)),
                              "PeacoQC" = peacoqc_kept,
                              "flowCut" = flowcut_kept,
                              "flowAI" = flowai_kept,
                              "flowClean" = flowclean_ketp)
    
    return(frame_cells)
}

ExtractChannels <- function(ff, remove_zeros){
    channels_no_scatter <- which(!is.na(pData(ff@parameters[,"desc"])))
    
    if (!remove_zeros){
        channels_no_scatter <- channels_no_scatter[grep("-A", 
                                                        ff@parameters@data[channels_no_scatter,"name"])]
    }
    
    if(length(grep(paste(c("FSC", "SSC"), 
                         collapse = "|"), 
                   ff@parameters@data[channels_no_scatter, "name"])) > 0){
        channels_no_scatter <- channels_no_scatter[-grep(paste(c("FSC", "SSC"), 
                                                               collapse = "|"), 
                                                         ff@parameters@data[channels_no_scatter, "name"])]}
    
    channels_scatter <- c(channels_no_scatter,
                          grep("FSC-A", ff@parameters@data[,"name"], fixed = TRUE),
                          grep("SSC-A", ff@parameters@data[,"name"], fixed = TRUE))
    
    channels_scatter <- unique(channels_scatter)
    
    channels_to_use_no_scatter <- colnames(ff)[channels_no_scatter]
    channels_to_use_scatter <- colnames(ff)[channels_scatter]
    
    if (any(c("Event_length","Beads", "Time", "Width") %in% pData(ff@parameters)[,"desc"])){
        
        names_channels <- pData(ff@parameters)[grep(paste(c("Event_length","Beads", "Time", "Width"), 
                                                          collapse = "|"), 
                                                    pData(ff@parameters)[,"desc"]),"name"]
        
        channels_to_use_no_scatter <- channels_to_use_no_scatter[!(channels_to_use_no_scatter %in% names_channels)]
        
        channels_to_use_scatter <- channels_to_use_scatter[!(channels_to_use_scatter %in% names_channels)]
        
    }
    
    return(list("no_scatter" = channels_to_use_no_scatter,
                "scatter" = channels_to_use_scatter))
    
}



Preprocessing <- function(ff, remove_margins = TRUE, comp_matrix, type = FALSE,
                          remove_zeros = FALSE){
    
    
    channels_to_use <- ExtractChannels(ff, remove_zeros)
    
    if (!is.null(comp_matrix)){
        comp_matrix <- ff@description[[comp_matrix]]}
    
    if (remove_margins == TRUE){
        ff <- RemoveMargins(ff, channels = channels_to_use$scatter)}
    if (!is.null(comp_matrix)){
        ff <- compensate(ff, comp_matrix)}
    
    if(!remove_zeros){
        if (length(channels_to_use$no_scatter) > 0){
            if (type == TRUE){
                ff <- transform(ff, estimateLogicle(ff, channels = channels_to_use$no_scatter,
                                                    type = "data"))
            }else{
                ff <- transform(ff, estimateLogicle(ff, channels = channels_to_use$no_scatter))
                
            }
        }
    } else {
        ff <- transform(ff, transformList(channels_to_use$no_scatter,
                                          arcsinhTransform(a = 0, b = 1/5, c = 0)))
    }
    return(ff)
}


flowrepo_quality_control <- function(ff, Figure_dir,
                                     remove_zeros = FALSE){
    channels_to_use <- ExtractChannels(ff, remove_zeros)
    
    
    tryout <- flowSet(ff)
    identifier(tryout[[1]]) <- basename(file)
    
    
    to_exclude_channels <- colnames(ff)[- grep(paste(channels_to_use$scatter, collapse = "|"), colnames(ff))]
    # Complete Default
    
    # flowai_full <- flow_auto_qc(tryout,
    #                             folder_results = file.path(Figure_dir, "FlowAI"),
    #                             output = 2, ChExcludeFS = "Original_ID",
    #                             ChExcludeFM = "Original_ID", mini_report = FALSE,
    #                             fcs_QC = FALSE, html_report = TRUE)
    # 
    # With figures
    # flowai_full <-  tryCatch({flow_auto_qc(tryout,
    #                                        folder_results = file.path(Figure_dir, "FlowAI_Selected_Channels"),
    #                                        output = 2, ChExcludeFS = to_exclude_channels,
    #                                        ChExcludeFM = to_exclude_channels, mini_report = FALSE,
    #                                        fcs_QC = FALSE, html_report = TRUE)},
    #                          error = function(e){return(NA)})

# Without figures    
    flowai_full <-  tryCatch({flow_auto_qc(tryout,
                                           folder_results = file.path(Figure_dir, "FlowAI_Selected_Channels"),
                                           output = 2, ChExcludeFS = to_exclude_channels,
                                           ChExcludeFM = to_exclude_channels, mini_report = FALSE,
                                           fcs_QC = FALSE, html_report = FALSE)},
                             error = function(e){return(NA)})
    # With figure
    # results_peacoqc <- PeacoQC(ff, 
    #                            channels = channels_to_use$scatter,
    #                            output_directory = 
    #                                Figure_dir,
    #                            name_directory = "PeacoQC_Results",
    #                            display_cells = 2000,
    #                            save_fcs = F,
    #                            plot = T,
    #                            remove_zeros = remove_zeros)
    # Without figure
    results_peacoqc <- PeacoQC(ff, 
                               channels = channels_to_use$scatter,
                               output_directory = 
                                   Figure_dir,
                               name_directory = "PeacoQC_Results",
                               display_cells = 2000,
                               save_fcs = F,
                               plot = F,
                               remove_zeros = remove_zeros)
    
    if (is.na(flowai_full)){
        frame_cells <- data.frame("idc" = seq_len(nrow(ff)),
                                  "PeacoQC" = results_peacoqc$GoodCells,
                                  "flowAI" = rep(NA, nrow(ff)))
    } else{
        
        frame_cells <- data.frame("idc" = seq_len(nrow(ff)),
                                  "PeacoQC" = results_peacoqc$GoodCells,
                                  "flowAI" = flowai_full@exprs[,"QCvector"] < 10000)
    }
    
    return(frame_cells)
    
}

flowrepo_overview_figure <- function(frame_cells, plot_channel, time_unit){
    
    frame_cells.melted <- melt(frame_cells,
                               measure.vars = c("PeacoQC",
                                                "flowAI"))
    if (nrow(ff) > 5000){
        sampled_cells <- sample(nrow(ff), 5000)}else{ 
            sampled_cells <- seq_len(nrow(ff))}
    
    channel <- plot_channel
    
    frame <- data.frame("idc" = sampled_cells,
                        "time" = ff@exprs[sampled_cells, "Time"],
                        "cells" = ff@exprs[sampled_cells,channel] )
    
    
    median_quantile_frame <-  CalculateMedianQuantile(ff, channel)
    if (nrow(ff) >=50000) {
        subset_timeplot <- sort(sample(seq_len(nrow(ff)), 50000))
    } else {
        subset_timeplot <- seq_len(nrow(ff))
    }
    
    
    h <- graphics::hist(flowCore::exprs(ff)[subset_timeplot, "Time"],
                        breaks=seq(min(flowCore::exprs(ff)[, "Time"]),
                                   max(flowCore::exprs(ff)[, "Time"]) +
                                       time_unit, by=time_unit),
                        plot=FALSE)
    
    h$events_1 <- sapply(h$mids, 
                         FUN = function(x){
                             which(abs(ff@exprs[,"Time"] - x) == 
                                       min(abs(ff@exprs[,"Time"]- x)))[1]})
    
    breaks <- sapply(h$breaks, 
                     FUN = function(x){
                         which(abs(ff@exprs[,"Time"] - x) ==
                                   min(abs(ff@exprs[,"Time"] - x)))[1]})
    
    
    scale_plot <- BuildManualScale(file)
    
    idcs <- findInterval(flowCore::exprs(ff)[subset_timeplot, "Time"], 
                         h$breaks)
    
    p_time <- ggplot() + theme_bw() +
        scale_plot +
        geom_point(aes(x = h$events_1, y = h$counts), size = 0.2) +
        geom_line(aes(x = h$events_1, y = h$counts)) +
        scale_x_continuous(limits = c(0,max(breaks))) +
        theme(panel.grid=element_blank(),
              axis.line  = element_line(size = 0.5),
              axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.y = element_line(size = 0.5),
              axis.ticks.x = element_blank(),
              panel.border = element_blank()) +
        xlab("") + ylab(paste0("Nr of cells per ",time_unit," seconds")) +
        theme(plot.title=element_text(size=10))
    
    p1 <- ggplot() +
        scale_plot +
        geom_point(data = frame, aes(x = idc, y = cells), size = 0.4) +
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.line.y = element_line(size = 0.5)) +
        labs(y = colnames(ff@exprs)[channel]) + 
        
        geom_line(data = median_quantile_frame,
                  aes(x = median_quantile_frame[,1], 
                      y = median_quantile_frame[,2]),
                  color = "red", size = 1)
    
    p2 <- HeatmapQCMassData(frame_cells.melted)
    
    legend <- get_legend(p2)
    legend_plot <- as_ggplot(legend)
    
    plot_without_legend <- ggarrange(p_time + theme(legend.position = "none"), 
                                     p1 + theme(legend.position = "none"), 
                                     p2 + theme(legend.position = "none"), 
                                     ncol = 1, 
                                     align = "v")
    
    final_plot <- ggarrange(plot_without_legend, legend_plot, 
                            nrow =1, widths = c(1,0.5))
    
    return(final_plot)
}


ExtractDatasetParameters <- function(dataset){
    # Dataset1 Flow
    
    if (dataset == "FR-FCM-Z2FV"){
        dataset_name <- "FR-FCM-Z2FV"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), 
                            pattern = ".fcs", 
                            full.names = TRUE)
        
        
        ff <- read.FCS(files[1])
        plot_channel <-  1
        time_unit <- 100
        comp_matrix <- "SPILL"
        remove_margins <- TRUE
        type <- FALSE
        remove_zeros <- FALSE
        
    }
    
    
    
    # Dataset2 # Containing histograms
    
    if (dataset == "FR-FCM-Z2XY"){
        dataset_name <- "FR-FCM-Z2XY"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".FCS", 
                            full.names = TRUE)
        channels <- c(1,3,5,6,8,9)
        
        
        plot_channel <-  1
        time_unit <- 100
        remove_zeros <- FALSE
        
    }
    
    # Dataset3 # Empty?? Flow
    
    
    if (dataset == "FR-FCM-Z2HL"){
        dataset_name <- "FR-FCM-Z2HL"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        files <- files[grep("cells", files)]
        
        plot_channel <-  7
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- NULL
        type <-  TRUE
        remove_zeros <- FALSE
        
    }
    
    # Dataset4 Spectral
    
    if (dataset == "FR-FCM-Z2QV"){
        dataset_name <- "FR-FCM-Z2QV"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        files <- files[grep("MC", files)]
        
        plot_channel <-  10
        time_unit <- 10000
        remove_margins <- FALSE
        comp_matrix <- NULL
        type <-  TRUE
        remove_zeros <- FALSE
        
    }
    
    
    # Dataset5 Spectral
    
    if (dataset == "FR-FCM-Z2EK"){
        dataset_name <- "FR-FCM-Z2EK"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        plot_channel <-  11
        time_unit <- 100000
        remove_margins <- FALSE
        comp_matrix <- NULL
        type <-  TRUE
        remove_zeros <- FALSE
    }
    
    
    # Dataset6 Spectral
    
    if (dataset == "FR-FCM-Z2EM"){
        dataset_name <- "FR-FCM-Z2EM"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        plot_channel <-  11
        time_unit <- 100000
        remove_margins <- FALSE
        comp_matrix <- NULL
        type <-  TRUE
        remove_zeros <- FALSE
        
    }
    
    # Dataset7 Spectral
    
    if (dataset =="FR-FCM-Z2EL" ){
        dataset_name <- "FR-FCM-Z2EL"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        plot_channel <-  11
        time_unit <- 100000
        remove_margins <- FALSE
        comp_matrix <- NULL
        type <-  TRUE
        remove_zeros <- FALSE
    }
    
    
    # Dataset8 Flow
    
    if (dataset == "FR-FCM-Z2KU"){
        dataset_name <- "FR-FCM-Z2KU"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        files <- files[grep("AC33", files)]
        
        plot_channel <-  4
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- "SPILL"
        type <- TRUE
        remove_zeros <- FALSE
        
    }
    
    # Dataset9 Mass
    
    if (dataset == "FR-FCM-Z2NR"){
        dataset_name <- "FR-FCM-Z2NR"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        plot_channel <-  28
        time_unit <- 100000
        remove_margins <- FALSE
        comp_matrix <- "SPILL"
        type <- TRUE
        remove_zeros <- TRUE
    }
    
    
    # Dataset10 Mass
    
    if (dataset == "FR-FCM-Z2NX"){
        dataset_name <- "FR-FCM-Z2NX"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        plot_channel <-  16
        time_unit <- 100000
        remove_margins <- FALSE
        comp_matrix <- "SPILL"
        type <- TRUE
        remove_zeros <- TRUE
    }
    
    
    # Dataset11 Flow
    
    if (dataset == "FR-FCM-Z2HR"){
        dataset_name <- "FR-FCM-Z2HR"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        plot_channel <-  c("CD42b", "CD2", "CD5")
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- "SPILL"
        type <- TRUE
        remove_zeros <- FALSE
        
    }
    
    # Dataset12 Flow
    
    if(dataset == "FR-FCM-Z28P"){
        dataset_name <- "FR-FCM-Z28P"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        plot_channel <-  3
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- "SPILL"
        type <- TRUE
        remove_zeros <- FALSE
        
    }
    
    # Dataset13 Flow
    
    if (dataset == "FR-FCM-Z28M"){
        dataset_name <- "FR-FCM-Z28M"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        files <- files[grep("Tyrodes", files)]
        
        plot_channel <-  3
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- "SPILL"
        type <-  TRUE
        remove_zeros <- FALSE
        
    }
    
    # Dataset14 Flow
    
    if (dataset == "FR-FCM-Z28L"){
        dataset_name <- "FR-FCM-Z28L"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        files <- files[grep("Tyrodes", files)]
        
        
        plot_channel <-  3
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- "SPILL"
        type <- TRUE
        remove_zeros <- FALSE
        
    }
    # Dataset15 Flow
    
    if (dataset == "FR-FCM-Z28J"){
        dataset_name <- "FR-FCM-Z28J"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        files <- files[grep("Tyrodes", files)]
        
        plot_channel <-  3
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- "SPILL"
        type <- TRUE
        remove_zeros <- FALSE
    }
    # Dataset16 Flow
    
    if (dataset == "FR-FCM-Z28W"){
        dataset_name <- "FR-FCM-Z28W"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        
        plot_channel <-  3
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- "SPILL"
        type <- TRUE
        remove_zeros <- FALSE
    }
    
    # Dataset17 Flow
    
    if (dataset == "FR-FCM-Z32U"){
        dataset_name <- "FR-FCM-Z32U"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        files <- files[grep("PBMC", files)]
        
        plot_channel <-  9
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- "SPILL"
        type <- TRUE
        remove_zeros <- FALSE
    }
    # Dataset18 Flow
    
    if (dataset == "FR-FCM-ZYRN"){
        dataset_name <- "FR-FCM-ZYRN"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        files <- files[grep("PBMC", files)]
        
        
        plot_channel <-  4
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- "SPILL"
        type <- TRUE
        remove_zeros <- FALSE
        
    }
    
    # Dataset19 Flow
    
    if (dataset == "FR-FCM-ZYRX"){
        dataset_name <- "FR-FCM-ZYRX"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        files <- files[grep("PBMC", files)]
        
        plot_channel <- 4
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- "SPILL"
        type <-  TRUE
        remove_zeros <- FALSE
    }
    
    # Dataset20 Flow # Only height parameters
    
    if (dataset == "FR-FCM-Z2FK"){
        dataset_name <- "FR-FCM-Z2FK"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        
        plot_channel <-  2
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- "SPILL"
        type <-  TRUE
        remove_zeros <-  FALSE
    }
    
    # Dataset21 Mass
    if (dataset == "FR-FCM-Z23C"){
        dataset_name <- "FR-FCM-Z23C"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        
        plot_channel <-  38
        time_unit <- 100000
        remove_margins <- FALSE
        comp_matrix <- NULL
        type <-  TRUE
        remove_zeros <- TRUE
    }
    # Dataset22 Flow
    
    if (dataset == "FR-FCM-Z2EY"){
        dataset_name <- "FR-FCM-Z2EY"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        files <- files[grep("Full staining", files)]
        
        plot_channel <-  36
        time_unit <- 20
        remove_margins <- TRUE
        comp_matrix <- NULL
        type <- TRUE
        remove_zeros <- FALSE
    }
    
    if (dataset == "FR-FCM-Z35Y"){
        # Dataset23 Flow
        dataset_name <- "FR-FCM-Z35Y"
        
        files <- list.files(paste0("Data_flowRepository/", dataset_name), pattern = ".fcs", 
                            full.names = TRUE)
        
        
        plot_channel <-  7
        time_unit <- 100
        remove_margins <- TRUE
        comp_matrix <- NULL
        type <-  TRUE
        remove_zeros <- FALSE
    }
    
    return(list("dataset_name" = dataset_name,
                "files" = files,
                "plot_channel" = plot_channel,
                "time_unit" = time_unit,
                "remove_margins" = remove_margins,
                "comp_matrix" = comp_matrix,
                "type" = type,
                "remove_zeros" = remove_zeros
    ))
    
    
}










