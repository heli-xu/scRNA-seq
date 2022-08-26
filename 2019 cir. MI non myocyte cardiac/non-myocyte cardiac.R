library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)
library(readr)

setwd("C:/Temple/scRNA/case study/non-myocyte cardiac")

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

save(NMCC_healthy, file = "results/Healthy_NMCC_cluster.rdata")

load("results/Healthy_NMCC_cluster.rdata")

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

save(TIP_MI,file = "results/TIP_MI_cluster.rdata")

load("results/TIP_MI_cluster.rdata")


##subsetting EC
TIP_MI_EC <- subset(TIP_MI, idents = c(1,9,10))

metadata_EC <- TIP_MI_EC@meta.data

Sham_cell <- metadata_EC %>% 
  rownames_to_column(var = "cell_id") %>% 
  filter(sample=="Sham") %>% 
  pull(cell_id)
  
MI_day7 <- metadata_EC %>% 
  rownames_to_column(var = "cell_id") %>% 
  filter(sample=="MI_day7") %>% 
  pull(cell_id)

#########################

DefaultAssay(TIP_MI) <- "RNA"

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


FeaturePlot(TIP_MI, 
            reduction = "tsne", 
            features = c("Itgam","Cd68", "Mertk","Itgax",
                         "Cd14","Fcgr1", "Ccr2",
                         "Timd4", "Cd163"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##2, 4, 7, 11, 14, 16, 17 11, 20, 21, 24 MP

##EC####
FeaturePlot(TIP_MI,
            reduction = "tsne",
            features = c("Pecam1","Icam2","Cdh5","Kdr"),
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

VlnPlot(TIP_MI, 
        features = c("Pecam1","Icam2","Cdh5"),
        pt.size = 0)

####fibroblasts####
FeaturePlot(TIP_MI,
            reduction = "tsne",
            features= c("Col1a1","Pdgfra","Tcf21","Postn"),
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

VlnPlot(TIP_MI, 
        features = c("Icam2","Col1a1","Pdgfra","Tcf21"),
        pt.size = 0)

#
###LEC markers####
VlnPlot(TIP_MI, 
        features = c("Pdpn","Flt4","Prox1","Ackr4","Msr1","Fcgr2b"),
        split.by = "sample",
        idents = c(1,9,10,15,18,23,24))

VlnPlot(TIP_MI, 
        features = c("Pecam1","Ptprc","Cd68","Adgre1","Pdpn"),
        pt.size = 0)

VlnPlot(TIP_MI, 
        features = c("Adgrg1","Adgrf5"),
        pt.size = 0)

VlnPlot(TIP_MI, 
        features = c("Icam2","Col1a1"),
        split.by = "sample",
        idents = 18)

#######explore ADGRs#########

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

###get conserved markers######
load("data/gene_annotations_mouse.rdata")

get_conserved <- function(cluster){
  FindConservedMarkers(TIP_MI,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    mutate(cluster_id = cluster)
}

conserved_markers <- map_dfr(c(1,9,10,11,15,18), get_conserved)

top10 <- conserved_markers %>% 
  mutate(avg_fc = (Sham_avg_logFC + MI_day3_avg_logFC+MI_day7_avg_logFC) /3) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc) %>% 
  ungroup()

####SCTransform (didn't use, took long)####

split_TIP <- SplitObject(TIP_MI, split.by = "sample")






##took~8min 
split_TIP <- split_TIP[c("Sham", "MI_day3","MI_day7")] %>%  
  map(~.x %>% 
        NormalizeData() %>% 
        SCTransform(verbose=FALSE))

integ_features <- SelectIntegrationFeatures(object.list = split_TIP, 
                                            nfeatures = 2000) 


# Prepare the SCT list object for integration
split_TIP <- PrepSCTIntegration(object.list = split_TIP, 
                               anchor.features = integ_features)

# CCA: Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_TIP, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features,
                                        verbose=FALSE)

# Integrate across conditions
TIP_integrated <- IntegrateData(anchorset = integ_anchors, 
                               normalization.method = "SCT",
                               verbose=FALSE)


TIP_integrated<- RunPCA(object = TIP_integrated,verbose = FALSE)


MI_integrated <- RunUMAP(MI_integrated,
                         dims = 1:30,
                         reduction = "pca")ElbowPlot(Object=NMCC_healthy, ndims = 40)

NMCC_healthy <- FindNeighbors(NMCC_healthy, dims = 1:30) %>% 
  FindClusters(resolution = 0.8) %>% 
  RunTSNE(reduction= "pca", dims =1:24)

DimPlot(NMCC_healthy, reduction = "tsne")

DefaultAssay(NMCC_healthy) <- "RNA"

##subsetting MP####

MI_MP <- subset(TIP_MI, idents = c(2, 4, 7, 11, 14, 16, 17, 11, 20, 21, 24)) %>%  
  NormalizeData() %>% 
  ScaleData() 

metadata_mp <- MI_MP@meta.data

Idents(MI_MP) <- metadata_mp$sample

levels(MI_MP) 

DotPlot(MI_MP,
        features = c( "Egfr"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

VlnPlot(MI_MP, features = "Egfr")

