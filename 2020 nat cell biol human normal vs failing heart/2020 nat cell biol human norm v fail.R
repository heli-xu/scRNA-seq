library(dplyr)
library(Seurat)
library(harmony)
library(tidyverse)
library(data.table)

library(ggplot2)
library(viridis)
library(readr)


####1. load raw data and cell info (metadata)####
##1st dataset: all normal##
counts <- fread("raw/healthy/GSE109816_normal_heart_umi_matrix.csv")

##2nd dataset: normal vs cHF, dHF##
counts2 <- fread("raw/norm_HF/human_heart_sc_umi.csv")
##much fewer normal cells, so had to combine datasets

#if merging objects w/ metadata, error about attribute not matching vector 
#if merging all counts, then create object, works, but a lot of NAs
counts_all <- counts2 %>% 
  left_join(counts, by="V1")
##keeping common features in both datasets

rownames(counts_all) <- counts_all$V1
#for matrix, rownames_to_column doesnt seem to work well

counts_all$V1 <- NULL

human_normal_HF <- CreateSeuratObject(counts_all, project = "human_normal_HF",
                                  min.cells =3, min.features  = 200)


####2. metadata wrangling####
metadata <- human_normal_HF@meta.data %>% 
  rownames_to_column("ID")

##it's a mess!!!!!###
#cell count discrepancy in cluster info; they don't add up to metadata count
#had to use barcode/sc info files
#we'll set condition first, celltype later

cell1 <- read_table("raw/healthy/GSE109816_normal_heart_cell_info.txt")
#it's a pretty messed up table, with column names all seperated wrong, 
#but for now we only ID and condition (all normal)
#read_table from readr can do some coerce, read.table less flexible==>error
cell1 <- cell1 %>% 
  select(ID, Type) %>% 
  mutate(group="Normal", 
         CellType = if_else(str_detect(Type, "_CM"), "CM", "NCM"))

cell2 <- read_table("raw/norm_HF/GSE121893_human_heart_sc_info_this is wrong. no CM.txt")
#yeah no CM, but we're not using that column
cell2 <- cell2 %>% 
  select(ID, Type) %>% 
  mutate(group = if_else(str_detect(Type, "HF_"), "HF", "Normal"),
         CellType = if_else(str_detect(Type, "_CM"), "CM", "NCM"))
  


# cell_info <- read.table("raw/healthy/GSE109816_normal_heart_cell_cluster_info.txt",
#                         header = T) %>%
#   select(ID, condition=Condition, ident= Cluster_ID, CellType) %>%
#   mutate(group="Normal")
# 
# 
# 
# cell_info2 <-  read.table("raw/norm_HF/GSE121893_all_heart_cell_cluster_info.txt",
#              header = T) %>%
#   select(ID, condition, ident) %>%
#  # filter(ID%in%cell_info$ID) %>%
#   # mutate(ID = paste0(ID,"_2")) %>%
#   mutate(group = if_else(str_detect(condition, "HF"), "HF", "Normal"),
#          CellType = case_when(str_detect(ident, "LA")~"CM",
#                               str_detect(ident, "LV")~"CM",
#                               str_detect(ident, "MP")~"MP",
#                               str_detect(ident, "FB")~"FB",
#                               str_detect(ident, "EC")~"EC",
#                               str_detect(ident, "SMC")~"SMC",
#                               str_detect(ident, "AV")~"CM"))
  


all_cell <- bind_rows(cell1, cell2)


metadata_full <- metadata %>% 
  right_join(all_cell, by= "ID") %>% 
  column_to_rownames("ID")

human_normal_HF@meta.data <- metadata_full

#Idents(human_normal) <- metadata_full$CellType

save(human_normal_HF, file = "data/human_normal_HF_metadata.rdata")
#no cell type, only CM or not


#cell_type <- read.table("raw/norm_HF/GSE121893_human_heart_sc_info.txt", header = T) %>% 
  #select(ID, CellType) %>% 
  #group_by(CellType) %>% 
  #count()
#this table is wrong, no CM included#



####Normalization without integration###no need to run 
human_HF_no_integ <- NormalizeData(human_HF) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData() %>% 
  RunPCA(npc = 30) 

DimPlot(human_HF_no_integ, reduction = "pca", group.by = "group")
##looks pretty separated
VlnPlot(human_HF_no_integ, features = "PC_1", group.by = "group")


####3. Harmony and Clustering####
load("data/human_normal_HF_metadata.rdata")

human_normal_HF_hmn <- NormalizeData(human_normal_HF, verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(npc=20, verbose= FALSE)

human_normal_HF_hmn <- human_normal_HF_hmn %>% 
  RunHarmony("group", plot_convergence = TRUE)

DimPlot(human_normal_HF_hmn, reduction = "harmony", group.by = "group")

VlnPlot(human_normal_HF_hmn, features = "harmony_1", group.by = "group")


###UMAP and clusters##

human_normal_HF_hmn <- human_normal_HF_hmn %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.6) 


metadata_hmn <- human_normal_HF_hmn@meta.data

Idents(human_normal_HF_hmn) <- metadata_hmn$CellType
#assign CM or not

DimPlot(human_normal_HF_hmn, reduction = "umap", label = T)

#double check CM markers
FeaturePlot(human_normal_HF_hmn,
            reduction = "umap", 
            features = c("MYH6", "RYR2", "TNNT2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

save(human_normal_HF_hmn, file = "data/human_normal_HF_hmn_cluster.rdata")


####4. Visualization######

######4.1 cell count####

##A little bit on cell count across different conditions
#used met.brewer for color pallet 
metadata_full %>% 
  ggplot(aes(x=group, fill=CellType)) +
  geom_bar() +
 # scale_fill_manual(values = met.brewer("Renoir",8, direction = -1, type = "discrete"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


######4.2 ADGRs of interest in CMs####
##subset CM###
load("data/human_normal_HF_hmn_cluster.rdata")

all_features <- rownames(human_normal_HF_hmn@assays$RNA@data)

human_CM <- subset(human_normal_HF_hmn, idents = "CM") %>% 
  NormalizeData() %>% 
  ScaleData(feature = all_features)

metadata_cm <- human_CM@meta.data

Idents(human_CM) <- metadata_cm$group

levels(human_CM) <- c("Normal", "HF")

DoHeatmap(human_CM, 
          features = c("MYH6","MYH7","ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                       "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))
#MYH6 MYH7 are LV and LA markers 

DotPlot(human_CM,
        features = c("MYH6","ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                     "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#adjust x axis label

DotPlot(human_CM,
        features = "ADGRF5",
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#without marker gene, dot will be more noticable

VlnPlot(human_CM, features = "ADGRF5")

CM_markers <- FindAllMarkers(human_CM, test.use = "MAST")

######4.3 ADGRs in ECs ####
load("data/human_normal_HF_hmn_cluster.rdata")

human_NCM <- subset(human_normal_HF_hmn, idents = "NCM")

metadata_NCM <- human_NCM@meta.data

Idents(human_NCM) <- metadata_NCM$seurat_clusters

#look at clusters by markers 
FeaturePlot(human_NCM,
            reduction = "umap", 
            features = c("PECAM1","VWF","FLT1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#1, 8, 11, 6 EC

FeaturePlot(human_NCM,
            reduction = "umap", 
            features = c("DCN", "CDH19"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##4 FB


FeaturePlot(human_NCM,
            reduction = "umap", 
            features = c("MYH11","MYL9","CASQ2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# SMC

FeaturePlot(human_NCM,
            reduction = "umap", 
            features = c("CSF1R", "CD14", "C1QC"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#Myeloid/MP: 10

FeaturePlot(human_NCM,
            reduction = "umap", 
            features = c("CD3E", "CD3D"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

human_EC <- subset(human_NCM, idents = c(1, 6, 8, 11))%>% 
  NormalizeData() %>% 
  ScaleData(feature = all_features)  #if not assigned, pull from upstream object


metadata_EC <- human_EC@meta.data

Idents(human_EC) <- metadata_EC$group

levels(human_EC) <- c("Normal", "HF")

DoHeatmap(human_EC, 
          features = c("PECAM1","ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                       "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))

VlnPlot(human_EC, features = "ADGRF5", split.by = "group")
VlnPlot(human_EC, features = "ADGRL4", split.by = "group")


DotPlot(human_EC,
        features = c("PECAM1",
          "ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                     "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

EC_markers <- FindAllMarkers(human_EC, test.use = "MAST")


####4.4 ADGRs in MP subset####
human_MP <- subset(human_NCM, idents = 10)%>% 
  NormalizeData() %>% 
  ScaleData(feature = all_features)


metadata_MP <- human_MP@meta.data

Idents(human_MP) <- metadata_MP$group

levels(human_MP) <- c("Normal", "HF")

DoHeatmap(human_MP, 
          features = c("CD14", #CD3D if for T cells
                       "ADGRF5", "ADGRG1", 
                       "ADGRE5", "ADGRE1", "ADGRE4P", #try both ADGRE4 and 4P
                       "ADGRA2", "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))

VlnPlot(human_MP, features = "ADGRF5", split.by = "group")

DotPlot(human_MP,
        features = c(#"CD14", #makes the rest of the dots very small
                     "ADGRF5", "ADGRG1", 
                     "ADGRE5", "ADGRE1", "ADGRE4P", #try both ADGRE4 and 4P
                     "ADGRA1", "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
