library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)

#########import

SHAM1.data <- Read10X(data.dir = "raw/Sham wk1")
SHAM1 <- CreateSeuratObject(counts=SHAM1.data,project = "TAC",
                              min.cells = 3,min.features = 200)

SHAM4.data <- Read10X(data.dir="raw/Sham wk4/")
SHAM4 <- CreateSeuratObject(counts=SHAM4.data,project = "TAC",
                            min.cells = 3,min.features = 200)

TAC1.data <- Read10X(data.dir="raw/TAC wk1/")
TAC1 <- CreateSeuratObject(counts=TAC1.data,project = "TAC",
                              min.cells = 3,min.features = 200)

TAC4.data <- Read10X(data.dir = "raw/TAC wk4/")
TAC4 <- CreateSeuratObject(counts=TAC4.data,project = "TAC",
                            min.cells = 3,min.features = 200)


###adding sample to metadata

SHAM1$sample <- sample("Sham1", size= ncol(SHAM1), replace = TRUE)
SHAM4$sample <- sample("Sham4", size= ncol(SHAM4), replace = TRUE)
TAC1$sample <- sample("TAC1", size= ncol(TAC1), replace = TRUE)
TAC4$sample <- sample("TAC4", size= ncol(TAC4), replace = TRUE)

head(x = SHAM1@meta.data, 5)

##merged datasets to a Seurat object
merged_TAC <- merge(x=SHAM1,
                   y=c(SHAM4, TAC1, TAC4),
                   add.cell.ids = c("SHAM1","SHAM4","TAC1","TAC4"))

metadata <- merged_TAC@meta.data

save(merged_TAC, file="data/merged_TAC.rdata")

load("data/merged_TAC.rdata")

metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


####Normalized without integration (group difference)
merged_TAC_no_integ <- NormalizeData(merged_TAC) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:25) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dim=1:25)

DimPlot(merged_TAC_no_integ, reduction = "umap", group.by = "sample")

save(merged_TAC_no_integ,file="data/merged_TAC_cluster.rdata")

load("data/merged_TAC_cluster.rdata")

####SCTransform with integration####

load("data/cell_cycle_mouse.rdata")

s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")


split_TAC <- SplitObject(merged_TAC, split.by = "sample")

####split by sample so remember to match case sensitive!
split_TAC <- split_TAC[c("Sham1","Sham4","TAC1","TAC4")] %>%  
  map(~.x %>% 
        NormalizeData() %>% 
        CellCycleScoring(g2m.features = g2m_genes, 
                         s.features = s_genes) %>% 
        SCTransform(verbose=FALSE))
##took a few min to run

save(split_TAC, file = "data/integrated_split_TAC.rdata")  ##47s
#saveRDS(split_TAC, file = "data/integrated_split_TAC2.rds") took the same time

integ_features <- SelectIntegrationFeatures(object.list = split_TAC, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_TAC <- PrepSCTIntegration(object.list = split_TAC, 
                               anchor.features = integ_features)

# CCA: Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_TAC, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features,
                                        verbose=FALSE)
#took 8m45s to run


# Integrate across conditions
TAC_integrated <- IntegrateData(anchorset = integ_anchors,
                               normalization.method = "SCT",
                               verbose=FALSE)
#took ~1.5min

###PCA and UMAP###
TAC_integrated <- RunPCA(object = TAC_integrated,verbose = FALSE) %>% 
  RunUMAP(dims = 1:30,
          reduction = "pca")

DimPlot(TAC_integrated, reduction = "umap", group.by = "sample")
#no cluster yet, just show overlapping among groups

##clusters###
##ElbowPlot(object = TAC_integrated, ndims = 40)##look for SD plataeu 

TAC_integrated <- FindNeighbors(object = TAC_integrated, 
                               dims = 1:20) %>%   ##1:20 based on paper
  FindClusters(resolution=0.8)

save(TAC_integrated, file = "data/clustered_integrated_TAC.rdata")

load("data/clustered_integrated_TAC.rdata")

DimPlot(TAC_integrated, reduction = "umap", label = TRUE, split.by = "sample")
#no segregation of clusters by samples

####Visualization of ADGRs in each cluster###
# Select the RNA counts slot to be the default assay
DefaultAssay(TAC_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
TAC_integrated <- NormalizeData(TAC_integrated, verbose = FALSE)

##briefly explore cell types, since it's not included in the object/metadata
#markers included in supplemental data of the paper
FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Cd14","Itgam","Cd68", "Ly6c2", 
                         "Ly6g", "Itgax"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##neutrophil 14

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Cd68","Fcgr1", "Adgre1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##macrophage 1, 2, 9, 10, 15, 17, 22, 23
VlnPlot(TAC_integrated, features = c("Cd68","Fcgr1", "Adgre1"), pt.size = 0 )

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Cd3d","Cd4","Cd8a", "Cd19"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##5,6,12,13,16 T cells
#0, 4,8, 18 B cells

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Gzma", "Klrb1c"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##3, 7 NK cells

##so we have half myeloid half lymphoid cells 

###Rename clusters####

TAC_integrated <- RenameIdents(object = TAC_integrated,
                               "0"="B",
                               "1"="macrophages",
                               "2"="macrophages",
                               "3"="NK",
                               "4"="B",
                               "5"="T",
                               "6"="T",
                               "7"="NK",
                               "9"="macrophages",
                               "10"="macrophages",
                               "12"="T",
                               "13"="B",
                               "14"="neutrophils",
                               "15"="macrophages",
                               "16"="T",
                               "17"="macrophages",
                               "18"="B",
                               "22"="macrophages",
                               "23"="macrophages")


FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Adgrg1","Adgrf5"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#very little

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Adgra1","Adgra2","Adgra3",
                         "Adgrb2","Adgrb3",
                         "Adgrc1", "Adgrc2","Adgrc3"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#very little
#not found: Adgra1, Adgrb3, Adgrc1, Adgrc2, Adgrc3
#Adgrb2 very little

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Adgrd1","Adgrd2",
                         "Adgre1","Adgre2", "Adgre3",
                         "Adgre4", "Adgre5"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# not found: Adgrd1, Adgrd3 (wrong--Adgrd2), Adgre2, Adgre3
#Adgre5 all clusters
#Adgre1: c2, C1, C15, C10some
#Adgre4: C10, C17

colors_in_use <- c("grey", "black", "orange", "red")
idents_in_use <- c(1,2,9,10,15,17,22,23)
VlnPlot(TAC_integrated, cols = colors_in_use, 
        idents= idents_in_use, features= "Adgre1", split.by = "sample")

VlnPlot(TAC_integrated, cols = colors_in_use,
        idents= idents_in_use, features= "Adgre4", split.by = "sample")

VlnPlot(TAC_integrated, cols = colors_in_use, pt.size=0, features= "Adgre5", split.by = "sample")

##to distinguish ccr2+/ccr2- macrophages
VlnPlot(TAC_integrated, pt.size = 0, idents = c(1,2,9,10,15,17,22,23), features = "Ccr2")
##seems like Adgre1(F4/80) and Adgre4 differentially expressed in different subsets of macrophages

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Adgrf1","Adgrf2", "Adgrf3",
                         "Adgrf4", "Adgrf5"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#very little
#not found: Adgrf1, Adgrf2, Adgrf4

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Adgrg1","Adgrg2", "Adgrg3",
                         "Adgrg4", "Adgrg5","Adgrg6","Adgrg7"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#not found: Adgrg2, Adgrg4,Adgrg7
#others not much
VlnPlot(TAC_integrated, idents= c(2,11,12), features= "Adgrg5", split.by = "sample")
VlnPlot(TAC_integrated, features= "Adgrg6")

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Adgrl1","Adgrl2", "Adgrl3",
                         "Adgrl4", "Adgrv1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#not much
VlnPlot(TAC_integrated, features= "Adgrl1")
VlnPlot(TAC_integrated, idents=c(1,2,5,6,12,13), features= "Adgrl1", split.by = "sample")


