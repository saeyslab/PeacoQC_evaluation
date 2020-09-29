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


files <- list.files("FlowCAP4/Flowframe_sampled/", pattern = ".fcs")

# Path to conversion file of FlowCAPIV files to svg names
rename <- read.table("FlowCAP4/SVG_files_manual_annotations/FlowCAP_svgNames.txt", 
                     header = F, 
                     sep = "")[,c(1,3)]
rename[,1] <- sprintf("%03d", rename[,1]) 


to_use_files <- rename[,1]


for(file in files[as.numeric(to_use_files)]){
 
     ff_new <- read.FCS(paste0("FlowCAP4/Flowframe_sampled/", file))
     
     expressionmatrix <- as_tibble(ff_new@exprs)
     
     colnames(expressionmatrix) <- paste0(colnames(expressionmatrix), 
                                          "_",pData(ff_new@parameters)[,2])
     
     expressionmatrix <- as_tibble(expressionmatrix)
     
     expressionmatrix <- as_tibble(cbind("rows" = as.numeric(rownames(expressionmatrix)), 
                                         expressionmatrix))
     
     
     expressionmatrix.long <- as_tibble(melt(expressionmatrix, 
                                             id.vars = c("Time_NA", "rows")))
     
     colnames(expressionmatrix.long) <- c("Time","event_ix", "channel_ix", "value")
     
     ncol <- 3
     nrow <- 2
     spacing_pt <- 30
     spacing <- unit(spacing_pt, "points")
     channel_width_pt <- 400
     channel_height_pt <- 300
     expand_x <-  0.05
     expand_y <-  0.05
     
     p <- ggplot(data = expressionmatrix.long, aes(y = value, x = event_ix)) +
         ggrastr::geom_point_rast() +
         facet_wrap( ~ channel_ix, scales = "free", ncol = ncol)+
         scale_x_continuous(expand = expand_scale(mult = expand_x,expand_x)) +
         scale_y_continuous(expand = expand_scale(mult = expand_y,expand_y)) +
         theme_void() +
         theme(panel.border = element_rect(fill = "#FFFFFF00", colour = "azure3"),
               strip.text.x = element_blank(),
               panel.spacing = spacing, legend.position = 'none')
     
     
     width <- ncol * channel_width_pt + (ncol-1) * spacing_pt
     height <- nrow * channel_height_pt + (nrow-1) * spacing_pt
     ggsave(paste0("FlowCAP4/Empty_SVG_plots/",
                   rename[as.numeric(substr(file, 1,3)),2],
                   ".svg"), 
            width = width/72.27, height = height/72.27)
     
     
     
}
