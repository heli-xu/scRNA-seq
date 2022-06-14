library(tidyverse)
library(Seurat)

library(R.utils)
##unzip data folder, but output is in current directory, not in raw/..
##or do it directly with 7-zip.;
#gunzip("raw/V1_Human_Heart_filtered_feature_bc_matrix.tar.gz", remove = FALSE)
#and then can use read10x to create object
#but this dataset from Genomics10x, not that clear what's what

#used dataset from the paper itself instead 
#downloaded from heartcellatlas.org, very big and comprehensive human heart data
#with multiple batch and different processing method, a lot of batch alignment using python
#downloaded data includes cell type info, BUT format in h5ad from Scanpy
#conversion from h5ad to Seurat object (not the count. including metadata and all)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

library(SeuratDisk)



Convert("raw/healthy_human_4chamber_map_unnormalized_V4.h5ad", dest = "h5seurat", overwrite = TRUE)

human_heart <- LoadH5Seurat("raw/healthy_human_4chamber_map_unnormalized_V4.h5seurat")
#less than 4min

save(human_heart, file = "data/2020human_heart_scRNA.rdata")


metadata2 <- human_heart@meta.data
##cannot use"$" for index

##already contains cell type data, plus cannot do normalization without batch alignment 
#and the integration isn't straightforward

Idents(human_heart) <- metadata2$cell_type
###IMPORTANT tip to add ident information, for visualization purposes
##looking at the object, there's reduction info (UMAP, PCA), so you can do DimPlot
##but no information of for cluster id


##For most visualization, set default back to RNA, although it might already be
##not sure if the object comes with normalization, just try vlnplot directly
#turns out no, if you do vlnplot without normalization
##expression levels is crazy high, and the dots are dashes within a cell type
DefaultAssay(human_heart) <- "RNA"

human_heart_normalized <- NormalizeData(human_heart, verbose = FALSE) 

save(human_heart_normalized, file="data/normalized_human_heart_scRNAseq.rdata")

#trying out one gene to look at plot
FeaturePlot(human_heart_normalized, 
            reduction = "umap", 
            features = c("ADGRF5"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

VlnPlot(human_heart_normalized, pt.size = 0, features = "ADGRF5")

feature_plot1 <- FeaturePlot(human_heart_normalized, 
            reduction = "umap", 
            features = c("ADGRA1","ADGRA2","ADGRA3",
                         "ADGRB2","ADGRB3",
                         "ADGRC1", "ADGRC2","ADGRC3",
                         "ADGRD1","ADGRD3",
                         "ADGRE1","ADGRE2", "ADGRE3"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##ADGRC1, 2,3 not found
##ADGRD3 not found 
##9 detectable genes in one figure for easier plot arrangemen

save(feature_plot1, file="plot/feature_plot1.rdata")

#vlnplot takes long, so we are not plotting unnecessary cluster (unassigned, doulets)
idents_allcell <- metadata2 %>%
  group_by(cell_type) %>%
  count() %>%
  filter(!cell_type %in% c("NotAssigned","doublets"))%>%
  pull(cell_type)

Vln_1 <- VlnPlot(human_heart_normalized, idents = idents_allcell,
                 pt.size = 0.01,
                 features = c("ADGRA2","ADGRA3","ADGRB3"))
#all present in CM, will be explored below
save(Vln_1, file = "plot/vln_a2a3b3.rdata")

feature_plot2 <- FeaturePlot(human_heart_normalized, 
            reduction = "umap", 
            features = c("ADGRE4P", "ADGRE5",
                         "ADGRF1","ADGRF2", "ADGRF3",
                         "ADGRF4", "ADGRF5"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# ADGRE4 not found

Vln_2 <- VlnPlot(human_heart_normalized, idents = idents_allcell,
                 pt.size = 0,
                 features = c("ADGRD1", "ADGRE1", "ADGRE2", "ADGRE5"))
##ADGRD1 present in CM, explored below in CM subsets
save(Vln_2, file= "plot/vln_d1e2e5")

VlnPlot(human_heart_normalized, idents=idents_allcell,
        pt.size =0,
        features = c("ADGRF3", "ADGRF5"))
#both have decent expression in CM, will be explored below

feature_plot3 <- FeaturePlot(human_heart_normalized, 
            reduction = "umap", 
            features = c("ADGRG1","ADGRG2", "ADGRG3",
                         "ADGRG4", "ADGRG5","ADGRG6"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

VlnPlot(human_heart_normalized, idents=idents_allcell,
        pt.size=0,
        features = c("ADGRG1","ADGRG2","ADGRG5", "ADGRG6"))
#all in CM, will be explored below

feature_plot4 <- FeaturePlot(human_heart_normalized, 
            reduction = "umap", 
            features = c("ADGRG7","ADGRL1","ADGRL2", 
                         "ADGRL3", "ADGRL4", "ADGRV1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

VlnPlot(human_heart_normalized, idents=idents_allcell,
        pt.size=0,
        features = c("ADGRL1", "ADGRL2", "ADGRL3", "ADGRL4"))
##all present in CM, ADGRL4 lowest, will be explored below 

##special interest##
VlnPlot(human_heart_normalized, idents = idents_allcell,
        features = c("ADGRF5","ADGRG1"))

##save feature plots together
#save(feature_plot1, feature_plot2, feature_plot3, feature_plot4,
     #file="plot/feature plots.rdata")
##took almost 5min

###subset vCM
human_heart_vCM <- subset(human_heart, idents = "Ventricular_Cardiomyocyte")

save(human_heart_vCM, file="data/human_heart_vCM_scRNAseq.rdata")

metadata_vCM <- human_heart_vCM@meta.data

Idents(human_heart_vCM) <- metadata_vCM$cell_states  

DefaultAssay(human_heart_vCM) <- "RNA"

DimPlot(human_heart_vCM, reduction = "umap")
##CM1 has largest number  

human_heart_vCM_normalized <- NormalizeData(human_heart_vCM)

save(human_heart_vCM_normalized, file="data/normalized_human_heart_vCM_scRNAseq.rdata")


##Featureplot doesn't do cluster very well, despite the labels, use vlnplot here 
VlnPlot(human_heart_vCM_normalized,
        pt.size = 0, features = c("ADGRA2","ADGRA3","ADGRB3"))

VlnPlot(human_heart_vCM_normalized, pt.size=0,
        features = c("ADGRD1","ADGRF3","ADGRF5"))

VlnPlot(human_heart_vCM_normalized, pt.size=0, 
        features = c("ADGRG1","ADGRG2","ADGRG5", "ADGRG6"))

VlnPlot(human_heart_vCM_normalized, pt.size=0,
        features = c("ADGRL1", "ADGRL2", "ADGRL3", "ADGRL4"))
