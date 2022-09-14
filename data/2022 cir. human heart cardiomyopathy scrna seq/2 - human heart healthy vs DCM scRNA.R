library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)
library(data.table)
library(readxl)
library(janitor)
library(stringr)
library(viridis)
library(harmony)

## 0. partition raw csv
#raw data file 23G, cannot generate seurat directly, 
#hence partition and generate new smaller files to save ram
#code in a different r file, only had to run once 

## 1. create seurat sub objects
#seurat objects generated and saved as rdata from partitioned data files 
#code in a different file, only had to run once

## 2. Import and merge
#need to adjust function, save as rds. or have to reassign to something else
M1 <- readRDS("processed data/seurat1.rds")
M2 <- readRDS("processed data/seurat2.rds")
M3 <- readRDS("processed data/seurat3.rds")
M4 <- readRDS("processed data/seurat4.rds")
M5 <- readRDS("processed data/seurat5.rds")
M6 <- readRDS("processed data/seurat6.rds")
M7 <- readRDS("processed data/seurat7.rds")
M8 <- readRDS("processed data/seurat8.rds")

human_dcm <- merge(M1, M2)
rm(M1, M2)
human_dcm <- merge(human_dcm, y= c(M3, M4))
rm(M3, M4)
human_dcm <- merge(human_dcm, y= c(M5, M6))
rm(M5, M6)
human_dcm <- merge(human_dcm, y= c(M7, M8))
rm(M7, M8)


## 3. metadata wrangling
#import metadata from supplemental
nuclei_data <- read_excel("raw/human dcm metadata.xlsx", sheet = 1) %>% 
  clean_names() %>% 
  select(samples, condition)

fresh_cell <- read_excel("raw/human dcm metadata.xlsx", sheet = 2) %>% 
  clean_names() %>% 
  rename(samples=sample) %>% 
  select(samples,condition)

metadata_condition <- bind_rows(nuclei_data, fresh_cell)

#metadata <- seurat1@meta.data
#rownames include single-cell id (sequence), and sample id in front of it
DCM_list <- metadata_condition %>% 
  filter(condition=="DCM") %>% 
  pull(samples)

#class(metadata) shows it's a dataframe but with rownames
#rownames usually get lost in dplyr cleaning, so make a column first
metadata_clean <- metadata %>% 
  rownames_to_column("cell_id") %>% 
  rowwise() %>% 
  mutate(condition = ifelse(str_detect(cell_id, paste(DCM_list, collapse = '|')), ##useful to detect multiple strings
                                    "DCM", "healthy")) %>% 
  column_to_ro

wnames("cell_id")  ##metadata format compliant

save(metadata_clean, file = "processed data/metadata_clean.rdata")

load("processed data/metadata_clean.rdata")


###4. insert metadata into merged obejct and normalization###
#metadata <- human_dcm@meta.data check columns

human_dcm@meta.data <- metadata_clean

human_dcm <- NormalizeData(human_dcm) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:15) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dim=1:20)

# human_dcm_hmn <- human_dcm %>% 
#   RunHarmony("condition", plot_convergence = TRUE)
# 
# human_dcm_hmn <- human_dcm_hmn %>% 
#   RunUMAP(reduction = "harmony", dims = 1:20) %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
#   FindClusters(resolution = 0.6) 
# 
# DimPlot(human_dcm_hmn, reduction = "umap", label = T)
#Run harmony change the looks of UMAP,
#but doesn't change the subset result (still has several clusters for one celltype)

save(human_dcm, file = "processed data/human_dcm_normalized_cluster.rdata")

###5. subset cardiomyocytes####
load("processed data/human_dcm_normalized_cluster.rdata")


metadata <- human_dcm@meta.data


##kinda not informative without cell types



DefaultAssay(human_dcm) <- "RNA"

FeaturePlot(human_dcm_hmn,
            reduction = "umap", 
            features = c("MYBPC3", "RYR2", "TNNT2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##2, 5 CM

VlnPlot(human_dcm, features = "RYR2", pt.size = 0)

FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("DCN", "CDH19"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##3, 4, 7, 11 FB

FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("PECAM1","VWF","FLT1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##0,8,9 EC



FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("RGS5","AGT","KCNJ8","PDGFRB"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##1 pericytes

FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("MYH11","MYL9","CASQ2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#12 smooth muscle

FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("CSF1R", "CD14", "C1QC", 
                         "FCGR1A", "CD68", "LY6G6C"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##6, 10, 13 myeloid
VlnPlot(human_dcm, features = c("ITGAM", "CD14","FCGR1A", "CD68"), pt.size = 0)

FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("NRXN1", "MPZ"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
VlnPlot(human_dcm, features = c("NRXN1", "MPZ"))
##17 neuronal cells

FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("CD3E", "CD2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("NKG7","KLRD1","GZMB"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##15, 14 T cells/ NK 

FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("CD79A","KIT"), 
            min.cutoff = 'q10', 
            label = TRUE)
VlnPlot(human_dcm, features = "KIT")
##not much B cells, 19 mast cells


FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("PLIN1","DGAT2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##20 adipocytes

FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("CCL21", "PDPN","PROX1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#16 lymphatics

FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("NRG1", "NRG3", "PCDH15"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##8 endocardial cells

FeaturePlot(human_dcm,
            reduction = "umap", 
            features = c("WWC1","PRG4","HAS1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

human_dcm <- RenameIdents(object = human_dcm,
                          "0"="EC",
                          "1"="Pericytes",
                              "2"="CM",
                              "3"="FB",
                              "4"="FB",
                              "5"="CM",
                              "6"="Myeloid",
                          "7" = "FB",
                          "8" = "EC",
                          "9" = "EC",
                          "10" = "Myeloid",
                          "11" ="FB",
                          "12" ="SMC",
                          "13" = "Myeloid",
                          "14" = "T cells/NK cells",
                          "15" = "T cells/NK cells",
                          "16" ="Lymphatics",
                          "17" = "Neurons",
                          "19" = "Mast cells",
                          "20" = "Adipocytes")

##umap cluster in the paper is neater, because harmonized 
#although the raw data says integrated, it doesn't mean 
#harmonized, just means it's combining the snRNA and scRNA
#CM doesn't get affected, because it's all from snRNA

cells_to_plot <- metadata %>% 
  rownames_to_column(var = "cell_id") %>% 
  filter(!CellType %in% c(18,21,22,23,24,25)) %>% 
  pull("cell_id")
##here I'm filtering out the scattered few cells without cell types

DimPlot(human_dcm, cells = cells_to_plot,
        reduction = "umap", label = T, label.size = 5)


####6. examining cell count and subsetting CM####
##look at cell count using metadata

human_dcm$CellType <- Idents(human_dcm)
##adding cell type to metadata for plotting

metadata$condition <- factor(metadata$condition, levels = c("healthy","DCM"))


metadata %>% 
  rownames_to_column(var = "cell_id") %>% 
  filter(cell_id %in% cells_to_plot) %>% ##only plotting cells with cell types
  ggplot(aes(x=condition, fill=CellType)) +
  geom_bar() +
  theme_minimal() +
  scale_fill_manual(values = DiscretePalette(24, palette = "stepped"))+
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


###subsetting  CM###
load("processed data/human_dcm_normalized_cluster.rdata")


human_dcm_CM <- subset(human_dcm, idents = "CM")

rm(human_dcm) #to save RAM

human_dcm_CM <- NormalizeData(human_dcm_CM) %>% 
  ScaleData()
##because we ran findvaribles (before scaling data) in the human_dcm, the default of 
#scaledata features are the variable features
##if there's no FindVariable before scaling, the default will be all features
#(but if we do want to include all features, see below)

metadata <- human_dcm_CM@meta.data

#now the idents is still 2, 5, we need to merge cluster, but split by condition
Idents(human_dcm_CM) <- metadata$condition
#check levels for plotting
levels(human_dcm_CM)

##check something from paper, to make sure subset is right 
DefaultAssay(human_dcm_CM) <- "RNA"
VlnPlot(human_dcm_CM, features = c("MYH6","ANKRD1","NPPA", "ADGRL3"), ncol=2, pt.size = 0)
#no need to split by condition


####7. exploring adgrs####


VlnPlot(human_dcm_CM, features = c("ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                                   "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"))

VlnPlot(human_dcm_CM, features = "ADGRF5")

DotPlot(human_dcm_CM,
        features = "ADGRF5",
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#if you only do genes expressed in similar% cells, the dots won't look too tiny,
#range will also change to only cover the selected genes

##scaling data with all features
#because otherwise some genes we want are not included in the variable features
all_features <- rownames(human_dcm_CM@assays$RNA@data)

human_dcm_CM_plot <- ScaleData(human_dcm_CM,
                               features = all_features)
###takes longer, and the objects ends up being much bigger (3G to 19G)

DoHeatmap(human_dcm_CM_plot, 
          features = c("MYH6","ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                       "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))



###8. subsetting EC####
load("processed data/human_dcm_normalized_cluster.rdata")

human_dcm_EC <- subset(human_dcm, idents = "EC")

rm(human_dcm) #to save RAM

human_dcm_EC <- NormalizeData(human_dcm_EC) %>% 
  ScaleData()

metadata_EC <- human_dcm_EC@meta.data

#we need to merge cluster, but split by condition
Idents(human_dcm_EC) <- metadata_EC$condition

levels(human_dcm_EC)

###visualization

DotPlot(human_dcm_EC,
        features = c("PECAM1",
                     "ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                     "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
        
VlnPlot(human_dcm_EC, features = "ADGRF5", pt.size = 0)  
VlnPlot(human_dcm_EC, features = "ADGRL4", pt.size = 0)
        
DoHeatmap(human_dcm_EC, 
          features = c("PECAM1","ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                       "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))        
        
####9. Subsetting MP####
load("processed data/human_dcm_normalized_cluster.rdata")

human_dcm_mye <- subset(human_dcm, idents = "Myeloid")

rm(human_dcm) #to save RAM

human_dcm_mye <- NormalizeData(human_dcm_mye) %>% 
  ScaleData()

metadata_mye <- human_dcm_mye@meta.data

#we need to merge cluster, but split by condition
Idents(human_dcm_mye) <- metadata_mye$condition

levels(human_dcm_mye) <- c("healthy", "DCM")

DotPlot(human_dcm_mye,
        features = c(#"CD14",
                     "ADGRF5", "ADGRG1", 
                     #"ADGRE5", 
                     "ADGRE1", "ADGRE4", "ADGRA1",
                     "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

