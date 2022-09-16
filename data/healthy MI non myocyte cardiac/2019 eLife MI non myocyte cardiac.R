library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)
library(readr)
library(harmony)
library(ggplot2)

####import data and clustering####

load("raw/sparce_matrix_TIP_ShamVsMI_days3_7.rdata")

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

save(TIP_MI,file = "clean/TIP_MI_cluster.rdata")
##cluster, no integration

load("clean/TIP_MI_cluster.rdata")

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

save(TIP_MI_hmn, file = "clean/TIP_MI_hmn_cluster.rdata")

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
##EC 2, 7, 9, 12, (15,) 18, 19

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
        features = c("Pdpn","Flt4","Prox1","Ackr4","Msr1","Fcgr2b", 
                     "Lyve1", "Fth1"),
       # split.by = "sample",
        idents = c(2, 7, 9, 12, 15, 18, 19))
##9, 15



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

###explore ADGRs in EC#########

ADGR_list <- c("Adgra1", "Adgra2","Adgra3", "Adgrb1","Adgrb2", "Adgrb3",
               "Adgrc3","Adgrc1", "Adgrc2", "Adgrd1", "Adgrd2", "Adgre1",
               "Adgre2", "Adgre3", "Adgre4", "Adgre5", "Adgrf1", "Adgrf2",
               "Adgrf3", "Adgrf4","Adgrf5","Adgrg1","Adgrg2","Adgrg3",
               "Adgrg4", "Adgrg5", "Adgrg6", "Adgrg7","Adgrl1", "Adgrl2",
               "Adgrl3", "Adgrl4", "Adgrv1")

VlnPlot(TIP_MI, idents=c(1,9,10,11,15,18),
        pt.size = 0,
        features= c("Ccl21a", "Mmrn1","Fgl2",
                    "Prss23","Thy1", "Igfbp5",
                    "Fth1","Lyve1","Prelp"),
        split.by = "sample")
#LEC markers


#####subsetting EC####
load("clean/TIP_MI_hmn_cluster.rdata")

TIP_MI_EC <- subset(TIP_MI_hmn, idents = c(2, 7, 9, 12, 15, 18, 19)) %>% 
  ScaleData() 

save(TIP_MI_EC, file = "clean/MI_EC_orig_ident.rdata")
##same object name as that with ident changed to sample (below), 
#so don't load them together

metadata_EC <- TIP_MI_EC@meta.data

Idents(TIP_MI_EC) <- metadata_EC$sample

levels(TIP_MI_EC) 

save(TIP_MI_EC, file = "clean/MI_EC_normalized.rdata")

DotPlot(TIP_MI_EC,
        features = c("Pecam1",
                     "Adgrf5", "Adgrg1", "Adgrl4", "Adgrl2", 
                     "Adgre5", "Adgre1", "Adgra2","Adgrd1"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

VlnPlot(TIP_MI_EC, features = c( "Adgrg1","Adgrf5"))

#####subset LEC####
load("clean/MI_EC_orig_ident.rdata")

VlnPlot(TIP_MI_EC, 
        features = c("Pdpn","Flt4","Prox1","Ackr4","Msr1","Fcgr2b", 
                     "Lyve1", "Fth1"),
        pt.size = 0.1)
        # split.by = "sample"

DotPlot(TIP_MI_EC,
        features = c("Pdpn","Flt4","Prox1","Ackr4","Msr1","Fcgr2b", 
                     "Lyve1"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


FeaturePlot(TIP_MI_EC, 
            reduction = "tsne", 
            features= ADGR_list, 
            sort.cell = TRUE, 
            min.cutoff = 'q10', 
            ncol = 4,
            label = TRUE)
#plotting in all EC to put it in perspective 

VlnPlot(TIP_MI_EC,
        features= c("Adgra2", "Adgra3", "Adgrb1", 
                    "Adgre1", "Adgre5", "Adgrf5", 
                    "Adgrg1", "Adgrg3", "Adgrl1", 
                    "Adgrl2", "Adgrl4"))


TIP_MI_LEC <- subset(TIP_MI_EC, idents = c(9, 15)) %>% 
  #NormalizeData() %>%  not needed after subset if normalization done before
  ScaleData() 
 

metadata_LEC <- TIP_MI_LEC@meta.data

Idents(TIP_MI_LEC) <- metadata_LEC$sample

levels(TIP_MI_LEC) 

save(TIP_MI_LEC, file = "clean/MI_LEC_normalized.rdata")



DotPlot(TIP_MI_LEC,
        features = c("Adgra2", "Adgra3", "Adgrb1", 
                       "Adgre1", "Adgre5", "Adgrf5", 
                       "Adgrg1", "Adgrg3", "Adgrl1", 
                       "Adgrl2", "Adgrl4"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

VlnPlot(TIP_MI_LEC, 
        features = c("Adgra2", "Adgra3","Adgrb1", 
                     "Adgre1", "Adgre5", "Adgrf5", 
                     "Adgrg1", "Adgrg3", "Adgrl1", 
                     "Adgrl2", "Adgrl4"),
        pt.size=0.1)

##if checking individual cluster###
VlnPlot(TIP_MI_EC, features =c("Adgra2", "Adgra3","Adgrb1", 
                                "Adgre1", "Adgre5", "Adgrf5", 
                                "Adgrg1", "Adgrg3", "Adgrl1", 
                                "Adgrl2", "Adgrl4"),
        idents = 9,
        group.by = "sample", #split.by doesn't show legend with multiple features
        pt.size=0.1)

#same as subsetting one single cluster and scale again

VlnPlot(TIP_MI_EC, features =c("Adgra2", "Adgra3","Adgrb1", 
                                 "Adgre1", "Adgre5", "Adgrf5", 
                                 "Adgrg1", "Adgrg3", "Adgrl1", 
                                 "Adgrl2", "Adgrl4"),
        idents = 15,
        group.by = "sample",
        pt.size=0.1)
####subsetting MP####

TIP_MI_MP <- subset(TIP_MI_hmn, idents = c(0, 4, 10, 11, 13, 14)) %>%  
  NormalizeData() %>% 
  ScaleData() 

metadata_mp <- TIP_MI_MP@meta.data

Idents(TIP_MI_MP) <- metadata_mp$sample

levels(TIP_MI_MP) 

save(TIP_MI_MP, file = "clean/MI_MP_normalized.rdata")

DotPlot(TIP_MI_MP,
        features = c(#"Cd14",
                     "Adgrf5", "Adgrg1", "Adgrl4", "Adgrl2", 
                     "Adgre5", "Adgre1", "Adgra2","Adgrd1"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

VlnPlot(MI_MP, features = "Egfr")

####saving all features as rdata####
all_features <- rownames(TIP_MI_hmn@assays$RNA@data)

save(all_features, file = "clean/all_features_MI.rdata")

