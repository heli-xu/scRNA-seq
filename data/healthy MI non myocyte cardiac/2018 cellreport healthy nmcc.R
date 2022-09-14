library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)
library(readr)



#NMCC.data <- read_csv("E-MTAB-6173.processed.1/full_count_matrix.txt")
#some error about stack usage
#use base R 
#Too big



load("data/sparce_matrix_fixed.rdata")

NMCC_healthy <- CreateSeuratObject(counts=sparce_matrix,project = "2017Skelly",
                                   min.cells = 3,min.features = 200)

metadata2 <- NMCC_healthy@meta.data

metadata2$cells <- rownames(metadata2)

metadata2 <- metadata2 %>% 
  mutate(sample = case_when(str_detect(cells,"-1")~"1",
                            str_detect(cells,"-2")~"2")) %>% 
  column_to_rownames(var = "cells")

NMCC_healthy@meta.data <- metadata2

metadata2 %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


NMCC_healthy <- NormalizeData(NMCC_healthy) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData() %>% 
  RunPCA() 

ElbowPlot(object=TIP_MI, ndims = 40)

NMCC_healthy <- FindNeighbors(NMCC_healthy, dims = 1:30) %>% 
  FindClusters(resolution = 1) %>% 
  RunTSNE(reduction= "pca", dims =1:24)

save(NMCC_healthy, file = "data/Healthy_NMCC_cluster.rdata")

load("data/Healthy_NMCC_cluster.rdata")

DimPlot(NMCC_healthy, reduction = "tsne")

NMCC_healthy.c4 <- subset(NMCC_healthy, idents = 4)
NMCC_healthy.c4 <- FindNeighbors(NMCC_healthy.c4, dims = 1:10)
NMCC_healthy.c4 <- FindClusters(NMCC_healthy.c4, resolution = 0.6) %>% 
  RunTSNE(reduction="pca",dims=1:10)

c4_markers <- FindAllMarkers(NMCC_healthy.c4,
                                  logfc.threshold = 0.25)

top10 <- c4_markers %>% 
  group_by(cluster) %>% 
  top_n(n=10,
        wt=avg_logFC) %>% 
  ungroup()

DimPlot(NMCC_healthy.c4, reduction = "tsne")



DefaultAssay(NMCC_healthy.c4) <- "RNA"


VlnPlot(NMCC_healthy.c4, 
        features = "Adgrg1")

VlnPlot(NMCC_healthy.c4, 
        features = "Adgrf5")


##Adhesion GPCR expression
FeaturePlot(NMCC_healthy.c4, 
            reduction = "tsne",
            features = c("Adgrg1","Adgrf5"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


##What is that cluster
FeaturePlot(NMCC_healthy,
            reduction = "tsne",
            features = c("Pecam1","Icam2","Cdh5"),
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
##looks like EC

FeaturePlot(NMCC_healthy,
            reduction = "tsne",
            features = c("Cd68","Il1b","H2-Aa","Cx3cr1"),
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)



