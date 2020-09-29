library(PeacoQC)
library(flowCore)


set.seed(1)

# Path to location of FlowCAP4 dataset
files <- list.files("../FlowRepositoryDatasets/FlowCAP4/", pattern = ".fcs")


for (file in files){
  ff <- read.FCS(paste0("../FlowRepositoryDatasets/FlowCAP4/", file))
  
  channels <- c(1,3,5,7,11,14)
  ff <- RemoveMargins(ff, channels)
  
  ff <- compensate(ff, ff@description$SPILL)
  
  if(ff@description$P3DISPLAY == "LOG"){
    
    ff <- transform(ff,transformList(c("SSC-A",colnames(ff@description$SPILL)), 
                                     logicleTransform()))
    
  } else{
    ff <- transform(ff,transformList(colnames(ff@description$SPILL), 
                                     logicleTransform()))
    
  }
  set.seed(1)
  
  sampled_cells <- sample(nrow(ff), 10000)
  
  ff_new <- ff[sort(sampled_cells),]
  
  write.FCS(ff, filename = paste0("../Inkscape/Flowframes_samples/", 
                                  sub(".fcs", "_Sampled.fcs", file)))
  
}
