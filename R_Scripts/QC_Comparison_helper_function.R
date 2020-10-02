# ---------------------------- Extract cells ----------------------------------

ExtractQCCells <- function(file, ff, directory){
    
    flowai_ff <- tryCatch({read.FCS(paste0(directory,"/FlowAI_full/",
                                           sub(".fcs", "_QC.fcs", file)))},
                          error = function(e){return(NA)})
    
    if (is(flowai_ff, "flowFrame")){
        flowai_kept <- flowai_ff@exprs[,"remove_from_all"] < 10000
        flowai_ff <- flowai_ff[which(flowai_ff@exprs[,"remove_from_all"] < 10000), ]
        cells_FlowAI <- ((nrow(ff) - nrow(flowai_ff))/nrow(ff))*100
        
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
                              "FlowAI" = flowai_kept,
                              "flowClean" = flowclean_ketp)
    
    return(frame_cells)
}


# ------------------------------- Calculate median and quantiles --------------

CalculateMedianQuantile <- function(ff, channels){
    breaks <- cut(seq_along(ff@exprs[,channels[1]]),
                  breaks = seq(0, nrow(ff@exprs)+1000, by = 1000),
                  labels = FALSE)
    mid_breaks <- seq(500, 
                      nrow(ff@exprs)+500, by = 1000)
    
    
    data_frame_list <- list()
    
    for (channel in channels){
    
    medians <- tapply(ff@exprs[,channel], breaks, median)
    q25 <- tapply(ff@exprs[,channel], breaks, quantile, 0.25)
    q75 <- tapply(ff@exprs[,channel], breaks, quantile, 0.75)
    
    
    data_frame_list[[colnames(ff@exprs)[channel]]] <- data.frame("idc" = mid_breaks,
                                        "median" = medians,
                                        "q25" = q25,
                                        "q75" = q75)
    }
    
    median_quantile_frame <- do.call(cbind, data_frame_list)
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
              panel.grid.minor = element_blank()) +
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
              panel.grid.minor = element_blank()) +
        scale_fill_manual(labels = c("FALSE" = "Removed by algorithm",
                                     "TRUE" = "Kept by algorithm",
                                     "NA" = "Algorithm did not work"),
                          values = c("TRUE" = "#27476E",
                                     "FALSE" = "#E26D5C",
                                     "NA" = "#DF2935")) 
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
        cells_FlowAI <- ((nrow(ff) - nrow(flowai_ff))/nrow(ff))*100
        
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
                              "FlowAI" = flowai_kept,
                              "flowClean" = flowclean_ketp)
    
    return(frame_cells)
}

