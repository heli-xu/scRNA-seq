#' Data Prep

library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)
library(data.table)



## Funciton to load CSV and export a Seurat
export_seurat = function(input, output){
  # input = inputs[8]; output = outputs[8]
  counts <- fread(input) %>% column_to_rownames("gene")
  seurat_tmp =  CreateSeuratObject(counts, project = "human DCM")
  saveRDS(seurat_tmp,file = output)
}

## get all inputs
inputs = list.files(path = "processed data/", pattern = "GSE183852_Integrated_Count", full.names = T)
inputs

## define all outputs
outputs =  paste0("processed data/seurat", 1:8 ,".rds")
outputs

#run one by one, to save RAM
export_seurat(inputs[1],outputs[1])
export_seurat(inputs[2],outputs[2])
export_seurat(inputs[3],outputs[3])
export_seurat(inputs[4],outputs[4])
export_seurat(inputs[5],outputs[5])
export_seurat(inputs[6],outputs[6])
export_seurat(inputs[7],outputs[7])
export_seurat(inputs[8],outputs[8])




## Loop!
#map2(inputs,outputs,~export_seurat(.x,.y))
