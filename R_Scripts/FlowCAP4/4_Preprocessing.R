library(PeacoQC)
library(flowCore)

files <- list.files("../FlowRepositoryDatasets/FlowCAP4/", pattern = ".fcs")

channels <- c(1,3,5,7,11,14)

for (file in files){
    ff <- read.FCS(file.path("../FlowRepositoryDatasets/FlowCAP4", file))
    
    ff <- RemoveMargins(ff, channels)
    
    ff <- compensate(ff, ff@description$SPILL)
    
    if(ff@description$P3DISPLAY == "LOG"){
        
        ff <- transform(ff,transformList(c("SSC-A",colnames(ff@description$SPILL)), 
                                         logicleTransform()))
        
    } else{
        ff <- transform(ff,transformList(colnames(ff@description$SPILL), 
                                         logicleTransform()))
        
    }  

    write.FCS(ff, paste0("FlowCAP4/Compensated_transformed_fcs/", 
                         sub(".fcs","_comp_trans.fcs",file)))
    
    }
