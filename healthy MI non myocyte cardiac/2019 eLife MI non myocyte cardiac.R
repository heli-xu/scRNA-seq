library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)
library(readr)
library(harmony)

####import data and clustering####

load("data/sparce_matrix_TIP_ShamVsMI_days3_7.rdata")

TIP_MI <- CreateSeuratObject(counts=sparce_matrix_TIP_ShamVsMI_days3_7,project = "2019Farbehi",
                             min.cells = 3,min.features = 200)

metadata <- TIP_MI@meta.data

metadata$cells <- rownames(metadata)

metadata <- metadata %>% 
  mutate(sample = case_when(str_detect(cells,"Sham")~"Sham",
                            str_detect(cells,"MI_day3")~"MI_day3",
                            str_detect(cells,"MI_day7")~"MI_day7")) %>% 
  column_to_rownames(var = "cells")

TIP_MI@meta.data <- metadata

metadata %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


TIP_MI <- NormalizeData(TIP_MI) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000, verbose = FALSE) %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(object=TIP_MI, ndims = 40)

TIP_MI <- FindNeighbors(TIP_MI, dims = 1:30) %>%
  FindClusters(resolution = 0.8) %>%
  RunTSNE(reduction= "pca", dims =1:24)

DimPlot(TIP_MI, reduction = "tsne", split.by = "sample")

save(TIP_MI,file = "data/TIP_MI_cluster.rdata")
##cluster, no integration

load("data/TIP_MI_cluster.rdata")

DimPlot(TIP_MI, reduction = "pca", group.by = "sample")
#not much separation by sample

###RunHarmony to integrate 
#not really changing clustering that much, but included here
TIP_MI_hmn <- TIP_MI %>% 
  RunHarmony("sample", plot_convergence = TRUE)

DimPlot(TIP_MI_hmn, reduction = "harmony", group.by = "sample")


TIP_MI_hmn <- TIP_MI_hmn %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(TIP_MI_hmn, reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.6)

DimPlot(TIP_MI_hmn, reduction = "tsne", split.by = "sample", label = T)

save(TIP_MI_hmn, file = "data/TIP_MI_hmn_cluster.rdata")

####id clusters by markers####

DefaultAssay(TIP_MI_hmn) <- "RNA"

FeaturePlot(TIP_MI,
            reduction = "tsne",
            features = c("Adgrg1","Adgra1"),
            split.by = "sample",
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

FeaturePlot(TIP_MI,
            reduction = "tsne",
            features = c("Col1a1","Postn","Cd68","Cd3d","Cd79a","Ccl5","Edr"),
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)


FeaturePlot(TIP_MI_hmn, 
            reduction = "tsne", 
            features = c("Itgam","Cd68", "Mertk","Itgax",
                         "Cd14","Fcgr1", "Ccr2",
                         "Timd4", "Cd163"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

VlnPlot(TIP_MI_hmn, 
        features = c("Itgam","Cd68","Cd14"),
        pt.size = 0)

##0, 4, 10, 11, 13, 14 MP


FeaturePlot(TIP_MI_hmn,
            reduction = "tsne",
            features = c("Pecam1","Icam2","Cdh5","Kdr"),
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)


VlnPlot(TIP_MI_hmn, 
        features = c("Pecam1","Icam2","Cdh5"),
        pt.size = 0)
##EC 2, 7, 9,  16, 18, 

FeaturePlot(TIP_MI_hmn,
            reduction = "tsne",
            features= c("Col1a1","Pdgfra","Tcf21","Postn"),
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

VlnPlot(TIP_MI_hmn, 
        features = c("Col1a1","Pdgfra","Tcf21"),
        pt.size = 0)
####fibroblasts 1,3, 5, 15, (16)

###LEC markers 
#not very clear
VlnPlot(TIP_MI_hmn, 
        features = c("Pdpn","Flt4","Prox1","Ackr4","Msr1","Fcgr2b"),
       # split.by = "sample",
        idents = c(2, 7, 9,  16, 18))

VlnPlot(TIP_MI_hmn, 
        features = c("Pecam1","Ptprc","Cd68","Adgre1","Pdpn"),
        pt.size = 0)

VlnPlot(TIP_MI_hmn, 
        features = c("Adgrg1","Adgrf5"),
        pt.size = 0)

VlnPlot(TIP_MI_hmn, 
        features = c("Icam2","Col1a1"),
        group.by = "sample",
        idents = 16)

#######explore ADGRs in EC#########
#need updating with TIP_MI_hmn and new idents

VlnPlot(TIP_MI, idents=c(1,9,10),
        features= "Adgra1", 
        split.by = "sample")

VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgra2", split.by = "sample")

VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgra3", split.by = "sample")

VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrb1", split.by = "sample")

VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrb2", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrb3", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrc3", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrc1", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrc2", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrd1", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrd2", split.by = "sample")

VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgre1", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgre2", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgre3", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgre4", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgre5", split.by = "sample")

VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrf1", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrf2", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrf3", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrf4", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrf5", split.by = "sample")

VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrg1", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrg2", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrg3", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrg4", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrg5", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrg6", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrg7", split.by = "sample")

VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrl1", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrl2", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrl3", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrl4", split.by = "sample")
VlnPlot(TIP_MI, idents=c(1,9,10),features= "Adgrv1", split.by = "sample")

VlnPlot(TIP_MI, idents=c(1,9,10,11,15,18),
        pt.size = 0,
        features= c("Ccl21a", "Mmrn1","Fgl2",
                    "Prss23","Thy1", "Igfbp5",
                    "Fth1","Lyve1","Prelp"),
        split.by = "sample")


VlnPlot(TIP_MI, idents=c(1,9,10,11,15,18,23,24),
        features= c("Adgre1","Adgre5","Adgrf5",
                    "Adgrg1","Adgrl2","Adgrl4"),
        pt.size = 0.2,
        split.by = "sample")

VlnPlot(TIP_MI, idents=c(1,9,10,11,15,18,23,24),
        features=c("Adgra2","Adgra3","Adgrg3","Cd68"),
        split.by = "sample")

VlnPlot(TIP_MI, idents=c(1,9,10,11,15,18),
        pt.size = 0,
        features=c("Adgra2","Adgra3","Adgre5","Adgrg3",
                   "Adgre1","Adgrf5",
                   "Adgrg1","Adgrl2","Adgrl4"),
        split.by = "sample")



##subsetting EC####
load("data/TIP_MI_hmn_cluster.rdata")

TIP_MI_EC <- subset(TIP_MI_hmn, idents = c(2, 7, 9,  16, 18)) %>% 
  NormalizeData() %>% 
  ScaleData() 

metadata_EC <- TIP_MI_EC@meta.data

Idents(TIP_MI_EC) <- metadata_EC$sample

levels(TIP_MI_EC) 

save(TIP_MI_EC, file = "data/MI_EC_normalized.rdata")

DotPlot(TIP_MI_EC,
        features = c( "Adgrg1","Adgrf5"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

VlnPlot(TIP_MI_EC, features = c( "Adgrg1","Adgrf5"))


##subsetting MP####

TIP_MI_MP <- subset(TIP_MI_hmn, idents = c(0, 4, 10, 11, 13, 14)) %>%  
  NormalizeData() %>% 
  ScaleData() 

metadata_mp <- TIP_MI_MP@meta.data

Idents(TIP_MI_MP) <- metadata_mp$sample

levels(TIP_MI_MP) 

save(TIP_MI_MP, file = "data/MI_MP_normalized.rdata")

DotPlot(MI_MP,
        features = c( "Egfr"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

VlnPlot(MI_MP, features = "Egfr")

####saving all features as rdata####
all_features <- rownames(TIP_MI_hmn@assays$RNA@data)

save(all_features, file = "data/all_features_MI.rdata")

