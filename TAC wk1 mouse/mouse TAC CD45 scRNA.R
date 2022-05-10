library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)

#########import 
##read10x file name has to be just barcodes/features/matrix.tsv.gz and nothing else

pool1.data <- Read10X(data.dir = "raw/GSE179276_pool1/")
pool1 <- CreateSeuratObject(counts=pool1.data$`Gene Expression`, project = "TAC wk1 CD45",
                            min.cells = 3, min.features = 200)

pool2.data <- Read10X(data.dir = "raw/GSE179276_pool2/")
pool2 <- CreateSeuratObject(counts=pool2.data$`Gene Expression`, project = "TAC wk1 CD45",
                            min.cells = 3, min.features = 200)

merged_TAC <- merge(x=pool1,
                    y=pool2,
                    add.cell.ids = c("pool1","pool2"))

#######Normalization without integration####
merged_TAC_no_integ <- NormalizeData(merged_TAC) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:25) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dim=1:25)

DimPlot(merged_TAC_no_integ, label = TRUE, reduction = "umap")

save(merged_TAC_no_integ, file = "data/merged_pools_TAC_no_integ.rdata")

###explore cell types###mostly interested in macrophages
DefaultAssay(merged_TAC_no_integ) <- "RNA"

FeaturePlot(merged_TAC_no_integ, 
            reduction = "umap", 
            features = c("Itgam","Cd68", "Ccr2","Timd4"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(merged_TAC_no_integ, 
            reduction = "umap", 
            features = c("Cd3d","Cd4","Cd8a", "Cd19"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
