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

load("data/merged_pools_TAC_no_integ.rdata")


###explore cell types###mostly interested in macrophages
DefaultAssay(merged_TAC_no_integ) <- "RNA"

FeaturePlot(merged_TAC_no_integ, 
            reduction = "umap", 
            features = c("Cd3g","Cd4","Cd8a", "Cd19","Cd79a"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#T cell B cell

FeaturePlot(merged_TAC_no_integ, 
            reduction = "umap", 
            features = c("Itgam","Cd68", "Mertk","Itgax"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#mono/macrophages

FeaturePlot(merged_TAC_no_integ, 
            reduction = "umap", 
            features = c("Cd14","Fcgr1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##more mono/macrophages

FeaturePlot(merged_TAC_no_integ, 
            reduction = "umap", 
            features = "Ly6g", 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#neutrophils

FeaturePlot(merged_TAC_no_integ, 
            reduction = "umap", 
            features = "Klrb1c", 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##NK cells

FeaturePlot(merged_TAC_no_integ, 
            reduction = "umap", 
            features = c("Timd4","Lyve1", "Cd163","Ccl24"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##resident macrophages

FeaturePlot(merged_TAC_no_integ, 
            reduction = "umap", 
            features = c("Ccr2", "H2-DMb1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##monocyte-derived macrophages (MoMF, recruited)

FeaturePlot(merged_TAC_no_integ, 
            reduction = "umap", 
            features = c("Ly6c2","Plac8"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##monocytes

##among monocyte/macrophages, adgre family has most expression
##cluster16 has low cd45 expression, but we'll throw it in there just incase 
VlnPlot(merged_TAC_no_integ, 
        features = "Adgre5", 
        pt.size = 0)

VlnPlot(merged_TAC_no_integ, 
        idents = c(0,1,2,3,5,9,11,13,16,17,23,25), 
        features = c("Adgre1","Adgre4"), 
        pt.size = 0)

##cross referencing resident macrophages/moMF
VlnPlot(merged_TAC_no_integ, 
        idents = c(0,1,2,3,5,9,11,13,16,17,23,25), 
        features = c("Ccr2"), pt.size = 0)

VlnPlot(merged_TAC_no_integ, 
        idents = c(0,1,2,3,5,9,11,13,17,23,25), 
        features = c("Timd4"), pt.size = 0)


##Adgrg6, Adgrl1 has some expression in macrophages, T cells respectively,
#quite low, violin plot not much showing up, 
##could be because of the lack of integration and filtering of double-hashtaged