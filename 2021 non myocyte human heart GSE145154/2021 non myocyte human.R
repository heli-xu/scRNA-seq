library(dplyr)
library(Seurat)
library(harmony)
library(tidyverse)
library(data.table)
library(purrr)
library(ggplot2)
library(viridis)
library(readr)
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

####1. import raw data and create objects####
##using PURRR::MAP
## list all files
folders = paste0("raw/",list.files("raw/"),"/")

data = folders %>% map(~Read10X(data.dir = .x))
S_objects = data %>% 
  map(~CreateSeuratObject(counts=.x, data.project= "human nonmyocyte",
                          min.cells = 3,min.features = 200))
#naming the object list
names(S_objects) <- list.files("raw/")

#adding sample to metadata
S_objects$`DCM-2-LVN`$sample <- sample("DCM", replace = TRUE)
S_objects$`DCM-2-LVP`$sample <- sample("DCM", replace = TRUE)
S_objects$`DCM-3-LVN`$sample <- sample("DCM", replace = T)
S_objects$`DCM-3-LVP`$sample <- sample("DCM", replace = T)
S_objects$`ICM-1-MIN`$sample <- sample("ICM", replace = T)
S_objects$`ICM-1-MIP`$sample <- sample("ICM", replace = T)
S_objects$`ICM-2-LVN`$sample <- sample("ICM", replace = T)
S_objects$`ICM-2-LVP`$sample <- sample("ICM", replace = T)
S_objects$`ICM-3-LVN`$sample <- sample("ICM", replace = T)
S_objects$`ICM-3-LVP`$sample <- sample("ICM", replace = T)
S_objects$`N-1-LVN`$sample <- sample("Normal", replace = T)
S_objects$`N-1-LVP`$sample <- sample("Normal", replace = T)

##merging objects

human_nmcc <- merge(x= S_objects[[1]], y= S_objects[2:12])
#note how grabbing first element is [[, but selecting multiple elements is []

save(human_nmcc, file="data/human_nmcc_unnormalized.rdata")

####2. Normalization without integration####
human_nmcc_no_integ <- NormalizeData(human_nmcc) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData() %>% 
  RunPCA(npc = 25) 

DimPlot(human_nmcc_no_integ, reduction = "pca", split.by = "sample")
##looks ok...
VlnPlot(human_nmcc_no_integ, features = "PC_1", group.by = "sample")

human_nmcc_cluster <- human_nmcc_no_integ %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dim=1:20)

save(human_nmcc_cluster, file="data/human_nmcc_normalized_clustered.rdata")

DimPlot(human_nmcc_cluster, reduction = "umap", label = T)


####Harmony and Clustering####no need##

# human_nmcc_hmn <- NormalizeData(human_nmcc, verbose = FALSE) %>% 
#   FindVariableFeatures(selection.method = "vst", 
#                        nfeatures = 2000, verbose = FALSE) %>% 
#   ScaleData(verbose = FALSE) %>% 
#   RunPCA(npc=20, verbose= FALSE)

human_nmcc_hmn <- human_nmcc_no_integ %>% 
  RunHarmony("sample", plot_convergence = TRUE)

DimPlot(human_nmcc_hmn, reduction = "harmony", split.by = "sample")

VlnPlot(human_nmcc_hmn, features = "harmony_1", group.by = "sample")
##looks pretty much the same with no_integration

###UMAP and clusters##

human_nmcc_hmn <- human_nmcc_hmn %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.6) 
###no need end here##




####3. Identification of cell types######

#look at clusters by markers 
FeaturePlot(human_nmcc_cluster,
            reduction = "umap", 
            features = c("PECAM1","VWF","FLT1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#5, 8, 11, 14, 18 EC

FeaturePlot(human_nmcc_cluster,
            reduction = "umap", 
            features = c("DCN", "CDH19"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##7, 13 FB


FeaturePlot(human_nmcc_cluster,
            reduction = "umap", 
            features = c("MYH11","MYL9","CASQ2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# 9, 10, 17* SMC

FeaturePlot(human_nmcc_cluster,
            reduction = "umap", 
            features = c("CSF1R", "CD14", "C1QC", "CD68",
                         "ITGAM", "FCGR1A"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#Myeloid/MP: 2,3,4,12,19
VlnPlot(human_nmcc_cluster, features = c("ITGAM", "CD14","FCGR1A", "CD68"), pt.size = 0)

FeaturePlot(human_nmcc_cluster,
            reduction = "umap", 
            features = c("CD3E", "CD3D",
                         "NKG7","KLRD1","GZMB"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##0,1,6, 21 T cell/NK cell

FeaturePlot(human_nmcc_cluster,
            reduction = "umap", 
            features = c("CD79A","KIT"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#16 B cell, 20 mast cells


####4. Subsetting EC####
all_features <- rownames(human_nmcc_cluster@assays$RNA@data)

human_nmcc_ec <- subset(human_nmcc_cluster, idents = c(5, 8, 11, 14, 18)) %>% 
  NormalizeData() %>% 
  ScaleData(feature = all_features)


metadata_EC <- human_nmcc_ec@meta.data

Idents(human_nmcc_ec) <- metadata_EC$sample

levels(human_nmcc_ec) <- c("Normal", "ICM", "DCM")

DoHeatmap(human_nmcc_ec, 
          features = c("PECAM1","ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                       "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))

VlnPlot(human_EC, features = "ADGRF5", split.by = "group")
VlnPlot(human_EC, features = "ADGRL4", split.by = "group")


DotPlot(human_nmcc_ec,
        features = c("PECAM1",
          "ADGRF5", "ADGRG1", 
          #"ADGRE5", "ADGRA1", "ADGRD1", "ADGRL1", "ADGRL3",
          "ADGRL4"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#EC_markers <- FindAllMarkers(human_EC, test.use = "MAST")


####5. Subsetting MP####
human_nmcc_MP <- subset(human_nmcc_cluster, idents = c(2,3,4,12,19))%>% 
  NormalizeData() %>% 
  ScaleData(feature = all_features)


metadata_MP <- human_nmcc_MP@meta.data

Idents(human_nmcc_MP) <- metadata_MP$sample

levels(human_nmcc_MP) <- c("Normal", "ICM", "DCM")

DoHeatmap(human_nmcc_MP, 
          features = c("CD14", #CD3D if for T cells
                       "ADGRF5", "ADGRG1", 
                       "ADGRE5", "ADGRE1", "ADGRE4P", #try both ADGRE4 and 4P
                       "ADGRA2", "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))

VlnPlot(human_MP, features = "ADGRF5", split.by = "group")

DotPlot(human_nmcc_MP,
        features = c(#"CD14", #makes the rest of the dots very small
                     "ADGRF5", "ADGRG1", 
                     "ADGRE5", 
                     "ADGRE1", "ADGRE4P", #try both ADGRE4 and 4P
                     "ADGRA1", "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
