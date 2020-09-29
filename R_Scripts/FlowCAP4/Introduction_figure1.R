library(flowCore)
library(ggplot2)
library(FlowSOM)
library(cowplot)


files <- list.files("FlowCAP4/Compensated_transformed_fcs/")

file <- files[147]

ff <- read.FCS(paste0("FlowCAP4/Compensated_transformed_fcs/", file))

channels <- c(5:14, 18, 20)

seed <- 1
set.seed(seed)


fsom <- FlowSOM(ff, colsToUse = channels, scale = FALSE)

labels <- cut(seq_len(nrow(ff)), 10, labels = FALSE)

test_colors <- grDevices::colorRampPalette(c("Lightblue", "blue"))(9)
test_colors <- c("red", test_colors)


subsampled_cells <- sort(sample(seq_len(nrow(ff)), 5000))

sampled_frame <- as.data.frame(ff@exprs[subsampled_cells,c("G710-A","R780-A", "Time")])

p_sample_1 <- ggplot() + geom_point(data = sampled_frame, aes(x = Time, y = `G710-A`), 
                                    color = test_colors[labels][subsampled_cells]) +
    theme_classic()

p_sample_2 <- ggplot() + geom_point(data = sampled_frame, aes(x = Time, y = `R780-A`), 
                                    color = test_colors[labels][subsampled_cells]) +
    theme_classic()

plot_time <- ggpubr::ggarrange(p_sample_1, p_sample_2)


angle <- pi/3
rotation_matrix <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 
                          nrow = 2)
fsom$MST$l <- t(apply(fsom$MST$l, 1, function(x) rotation_matrix %*% x))  

fsom_mirrored <- fsom
fsom_mirrored$MST$l[,2] <- -fsom_mirrored$MST$l[,2]

plot <- PlotPies(fsom_mirrored, cellTypes = labels, maxNodeSize = 0.2, 
                 colorPalette = test_colors)


plot_scatters <- Plot2DScatters(fsom, channelpairs = list(c("Time", "G710-A"), 
                                                          c("Time", "R780-A")), 
                                cluster = list(c(81,82,91,92, 93), 
                                               c(87,88,90,99,100)), plotFile = NULL)



arranged_scatters <- ggpubr::ggarrange(plotlist =  plot_scatters, 
                                       ncol = 2, nrow = 2)



full_plot <- ggpubr::ggarrange(ggpubr::ggarrange(plot_time, plot + 
                                                     theme(legend.position = "none"), 
                                                 nrow = 2),
                               arranged_scatters,
                               ncol = 2) 

ggsave(plot = full_plot, file = "FlowCAP4/Introduction_figure.png", width = 12, height = 10)
