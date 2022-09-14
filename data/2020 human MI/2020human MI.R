library(tidyverse)
library(Seurat)


if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")



library(SeuratDisk)
library(Matrix)
library(viridis)

#downloaded data includes cell type info, BUT format in h5ad from Scanpy
#conversion from h5ad to Seurat object (not the count. including metadata and all)


Convert("raw/snRNA-seq-submission.h5ad", dest = "h5seurat", overwrite = TRUE)

hfile<- Connect("raw/snRNA-seq-submission.h5seurat")

x <- hfile[["assays/RNA/data"]]
sp <- sparseMatrix(i=x[["indices"]][]+1,p=x[["indptr"]][],x=x[["data"]][]  )
dim(sp)
h5attr(x = x, which = "dims")

human_MI <- LoadH5Seurat("raw/snRNA-seq-submission.h5seurat") ## error

#Error in sparseMatrix????????????
#had to modify source code of LoadH5Seurat function (in separate R scripts)

human_MI_diy = object; save(human_MI_diy, file = 'human_MI_diy.rdata')
#object returned from "DIY count and LoadH5seurat modified.r' 

metadata <- human_MI_diy@meta.data
##cannot use"$" for index

##already contains cell type data, 

Idents(human_MI_diy) <- metadata$cell_type_original
###IMPORTANT tip to add ident information, for visualization purposes
##looking at the object, there's reduction info (UMAP, PCA), so you can do DimPlot

DimPlot(human_MI_diy,
        reduction = "umap",
        label = T,
        split.by = "major_labl")
#cell typess and cell count across different conditions are quite comparable


###subset CM, scale data within CM with all features
all_features <- rownames(human_MI_diy@assays$RNA@data)

human_MI_diy_CM <- subset(human_MI_diy, idents = "Cardiomyocyte")

human_MI_diy_CM <- NormalizeData(human_MI_diy_CM) %>% 
  ScaleData(feature = all_features)
#increase object size dramatically

#save(human_MI_diy_CM, file="data/human_MI_CM_non-normalized.rdata")

metadata_CM <- human_MI_diy_CM@meta.data

Idents(human_MI_diy_CM) <- metadata_CM$major_labl  

levels(human_MI_diy_CM) <- c("CTRL", "IZ","FZ", "BZ", "RZ")



DoHeatmap(human_MI_diy_CM, 
          features = c("MYH6","MYH7","ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                       "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))

DotPlot(human_MI_diy_CM,
        features = c("MYH6","ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                     "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

VlnPlot(human_MI_diy_CM, features = "ADGRF5")


###subsetting EC 
human_MI_EC <- subset(human_MI_diy, idents = "Endothelial")%>% 
  NormalizeData() %>% 
  ScaleData(feature = all_features)

metadata_EC <- human_MI_EC@meta.data

Idents(human_MI_EC) <- metadata_EC$major_labl  

levels(human_MI_EC) <- c("CTRL", "IZ","FZ", "BZ", "RZ")

DoHeatmap(human_MI_EC, 
          features = c("PECAM1","ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                       "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))

DotPlot(human_MI_EC,
        features = c("PECAM1","ADGRF5", "ADGRG1", "ADGRE5", "ADGRA1",
                     "ADGRD1", "ADGRL1", "ADGRL3","ADGRL4"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

VlnPlot(human_MI_EC, features = "ADGRF5", pt.size = 0)
