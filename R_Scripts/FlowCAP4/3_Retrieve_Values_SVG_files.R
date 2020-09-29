library(flowCore)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)
library(reshape2)
library(tidyverse)
library(xml2)
library(ggrastr) # to quickly plot 
library(ggplot2)
library(svglite)
library(pdist)

# Path to conversion file of FlowCAPIV files to svg names
rename <- read.table("FlowCAP4/SVG_files_manual_annotations/FlowCAP_svgNames.txt", 
                     header = F, 
                     sep = "")[,c(1,3)]
rename[,1] <- sprintf("%03d", rename[,1]) 

group <- 1
for (group in c(1:5) ){
    
    
    # Path to svg files
    files <- list.files("FlowCAP4/SVG_files_manual_annotations", 
                        full.names = T, 
                        pattern =  paste0("Group",group,".*svg"), 
                        recursive = TRUE)
    
    for(file in 1:11){
        
        manual_file <- paste0("Group",group,"_",file,  ".svg") 
        
        manual_files <- files[grep(manual_file, files)]
        
        
        # Path to sampled data that was used for creation of svg files
        ff_new <- read.FCS(paste0("FlowCAP4/Flowframe_sampled/", 
                                  rename[rename[,2] == sub(".svg", "", manual_file),1], "_Sampled.fcs"))
        
        colnames(ff_new@exprs) <- c(1:ncol(ff_new@exprs))
        
        expressionmatrix <- as.matrix(ff_new@exprs[,1:6])
        
        expressionmatrix <- cbind("rows" = as.numeric(rownames(expressionmatrix)), 
                                  expressionmatrix)
        
        
        expressionmatrix.long <- melt(expressionmatrix, 
                                      id.vars = "rows", 
                                      factorsAsStrings = F)
        
        colnames(expressionmatrix.long) <- c("event_ix", "channel_ix", "value")
        
        
        # Use exact same specifics of the ones used for creation of svg files
        ncol <- 3
        nrow <- 2
        spacing_pt <- 30
        spacing <- unit(spacing_pt, "points")
        channel_width_pt <- 400
        channel_height_pt <- 300
        expand_x <-  0.05
        expand_y <-  0.05
        n_channels <- 6
        n_events <- 10000
        
        channels <- tibble(
            channel_ix = 1:n_channels,
            n_events = n_events) %>%
            mutate(
                col = ((channel_ix - 1) %% ncol) + 1,
                row = ceiling(channel_ix / ncol),
                x = (col - 1) * (channel_width_pt + spacing_pt) + channel_width_pt * expand_x,
                y = (row - 1) * (channel_height_pt + spacing_pt)+ channel_height_pt * expand_y,
                width = channel_width_pt - 2 * channel_width_pt * expand_x,
                height = channel_height_pt - 2 * channel_height_pt * expand_y
            ) %>%
            mutate(
                y_mean = y + height / 2,
                x_mean = x + width / 2
            )
        
        
        end_data <- list()
        
        for (manual_annotated_file in manual_files){
            
            annotator <- strsplit(manual_annotated_file, "/Group._")[[1]][2]
            
            
            # get drawn rectangles
            svg <- read_xml(manual_annotated_file)
            rects <- xml_find_all(svg, ".//svg:rect")
            drawn_rects_filter <- rects %>% 
                xml_attr("style") %>% 
                stringr::str_detect("fill:#f90521")
            drawn_rects_filter[is.na(drawn_rects_filter)] <- FALSE
            drawn_rect_elements <- rects[drawn_rects_filter]
            
            
            # put in a nice tibble with x, y, width and height
            extract_rect_coords <- function(rect) {
                tibble(
                    width = as.numeric(xml_attr(rect, "width")),
                    height = as.numeric(xml_attr(rect, "height")),
                    x = as.numeric(xml_attr(rect, "x")),
                    y = as.numeric(xml_attr(rect, "y"))
                )
            }
            drawn_rects <- drawn_rect_elements %>% map_dfr(extract_rect_coords)
            
            if (nrow(drawn_rects) > 0){
                drawn_rects <- drawn_rects %>% 
                    mutate(
                        y_mean = y + height / 2,
                        x_mean = x + width / 2
                    ) %>% 
                    mutate(rect_ix = row_number())
                
                # add affected channel
                map_to_channel <- function(x, y, x_ref, y_ref) {
                    check <- matrix(c(x, y), byrow=FALSE, ncol = 2)
                    ref <- matrix(c(x_ref, y_ref), byrow=FALSE, ncol = 2)
                    
                    as.matrix(pdist::pdist(check, ref)) %>% 
                        apply(1, which.min)
                }
                
                drawn_rects <- drawn_rects %>% 
                    mutate(channel_ix = map_to_channel(x_mean, 
                                                       y_mean, 
                                                       channels$x_mean, 
                                                       channels$y_mean))
                
                # map to channel information and extract selected events
                drawn_rects_mapped <- drawn_rects %>% 
                    left_join(channels, "channel_ix", suffix = c("", "_channel")) %>% 
                    mutate(
                        x_local = x - x_channel,
                        percentage_start = x_local / width_channel,
                        percentage_end = (x_local + width) / width_channel,
                        event_ix_start = floor(percentage_start * n_events) %>% pmax(0),
                        event_ix_end = floor(percentage_end * n_events) %>% pmin(n_events)
                    )
                
                # create a mask of every affected event
                data_masked <- drawn_rects_mapped %>% 
                    select(channel_ix, rect_ix, event_ix_start, event_ix_end) %>% 
                    mutate(event_ix = map2(event_ix_start, event_ix_end, seq)) %>% 
                    unnest(c(event_ix)) %>% 
                    mutate(masked = TRUE)
                
                data_result <- expressionmatrix.long %>% 
                    full_join(data_masked, c("channel_ix", "event_ix")) %>% 
                    group_by(channel_ix, event_ix) %>% 
                    summarise(
                        masked = any(!is.na(rect_ix)),
                        rect_ix = list(rect_ix[!is.na(rect_ix)])
                    ) %>% 
                    ungroup() %>% 
                    left_join(expressionmatrix.long, c("channel_ix", "event_ix")) %>% 
                    filter(event_ix > 0)
                
                data_result_test <- data_result %>%
                    group_by(event_ix) %>% 
                    summarise(!!annotator := any(masked))
                
                
            } else {
                
                data_result_test <- as_tibble(data.frame("event_ix" = c(1:10000), 
                                                         annotator = rep(FALSE, 10000)))
                colnames(data_result_test) <-  c("event_ix", annotator)
                
             
            }
            end_data[[annotator]] <- data_result_test
            
            
            
        }
        full_data <- Reduce(full_join, end_data)
        
        
        # Path to folder where results can be stored
        saveRDS(full_data, 
                file = paste0("FlowCAP4/Full_annotated_tables/group", 
                              group, "_", file,".Rds"))
        
    }
}


