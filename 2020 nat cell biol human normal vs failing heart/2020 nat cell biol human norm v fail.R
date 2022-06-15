library(dplyr)
library(Seurat)
library(harmony)
library(tidyverse)
library(data.table)
library(MetBrewer)
library(ggplot2)
library(viridis)

###load raw data and cell info (metadata)##
###1. normal####
counts <- fread("raw/healthy/GSE109816_normal_heart_umi_matrix.csv")

rownames(counts) <- counts$V1
#for matrix, rownames_to_column doesnt seem to work well

counts$V1 <- NULL
##800m, doable, but stuck at next step

human_normal <- CreateSeuratObject(counts, project = "human normal",
                                  min.cells =3, min.features  = 200)
save(TAC_CM_NMCC, file = "data/TAC_CM_NMCC.rdata")

##too big, although the object turns out to be 400m


###2. normal vs cHF, dCM####
counts2 <- fread("raw/norm_HF/human_heart_sc_umi.csv")

rownames(counts2) <- counts2$V1

counts2$V1 <- NULL

human_HF <- CreateSeuratObject(counts2, project = "human HF",
                                   min.cells =3, min.features  = 200)
##i think it meant to be combined with that other healthy dataset


cell_info <- read.table("raw/norm_HF/GSE121893_all_heart_cell_cluster_info.txt", header = T)

#cell_type <- read.table("raw/norm_HF/GSE121893_human_heart_sc_info.txt", header = T) %>% 
  #select(ID, CellType) %>% 
  #group_by(CellType) %>% 
  #count()
#this table is wrong, no CM included#


###metadata wrangling####

metadata <- human_HF@meta.data %>% 
  rownames_to_column(var = "cell_id")

cell_info <- cell_info %>% 
  rename(cell_id = ID)

metadata_full <- metadata %>% 
  left_join(cell_info, by = "cell_id") %>% 
  mutate(group = if_else(str_detect(condition, "HF"), "HF", "Normal"),
         CellType = case_when(str_detect(ident, "LA")~"CM",
                              str_detect(ident, "LV")~"CM",
                              str_detect(ident, "MP")~"MP",
                              str_detect(ident, "FB")~"FB",
                              str_detect(ident, "EC")~"EC",
                              str_detect(ident, "SMC")~"SMC",
                              str_detect(ident, "AV")~"CM")) %>% 
  #filter(is.na(CellType))   #check for missing celltype
  column_to_rownames("cell_id")
##keep the rownames as cell_id

##metadata_full$condition <- factor(metadata_full$condition, levels = c("0w", "2w","5w","8w","11w"))
##good extra thing to do for future visualization

human_HF@meta.data <- metadata_full

Idents(human_HF) <- metadata_full$ident

save(human_HF, file = "data/human_HF_metadata.rdata")
##conditions no levels, cell type added


####Normalization without integration###no need to run 
human_HF_no_integ <- NormalizeData(human_HF) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData() %>% 
  RunPCA(npc = 30) 

DimPlot(human_HF_no_integ, reduction = "pca", group.by = "group")
##looks pretty separated
VlnPlot(human_HF_no_integ, features = "PC_1", group.by = "group")


####Harmony ####
human_HF_harmony <- NormalizeData(human_HF, verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(npc=20, verbose= FALSE)

human_HF_harmony <- human_HF_harmony %>% 
  RunHarmony("group", plot_convergence = TRUE)

DimPlot(human_HF_harmony, reduction = "harmony", group.by = "group")

VlnPlot(human_HF_harmony, features = "harmony_1", group.by = "group")


####UMAP and clusters####

human_HF_harmony <- human_HF_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.6) 


DimPlot(human_HF_harmony, reduction = "umap", label = T, split.by = "group")

Idents(human_HF_harmony) <- metadata_full$CellType

##run normalization and umap for visualization
TAC_integrated <- RunPCA(object = TAC_integrated,verbose = FALSE) %>% 
  RunUMAP(dims = 1:20,
          reduction = "pca")


DimPlot(TAC_integrated, reduction = "umap", group.by = "condition")

save(TAC_integrated, file = "data/TAC_integrarted.rdata")

####################################################################
###Visualization of clusters###
load("data/TAC_integrarted.rdata")

metadata <- TAC_integrated@meta.data 

##rearrange condition levels, for easier visualization
metadata$condition <- factor(metadata$condition, levels = c("0w", "2w","5w","8w","11w"))

#insert back to object
TAC_integrated@meta.data <- metadata

##label whatever object you're plotting
Idents(TAC_integrated) <- metadata$CellType


DimPlot(TAC_integrated, reduction = "umap", label = T)

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
##cell count 2000-3000, majority CM


###Explore ADGRs###

DefaultAssay(TAC_integrated) <- "RNA"

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Gpr123","Gpr124","Gpr125",
                         "Bai1","Bai2",
                         "Cselr1", "Cselr2","Cselr3"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##somehow this paper doesn't use new names, but the old ones, 
#so none of the adgrs showed up, had to search the gprs, and some other names
#Gpr124 mainly in FB

VlnPlot(TAC_integrated, pt.size = 0.1, idents = idents_to_use,
        features = "Gpr124", split.by = "condition")
#not much 

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Gpr133","Gpr144",
                         "Emr1","Emr2", "Emr3",
                         "Emr4", "Cd97"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#assign a pallet for violin plot
color_timepoint <- met.brewer("Homer1", 5, direction = -1, type = "discrete")
idents_to_use <- c("CM","MP","T","FB","EC","GN")

VlnPlot(TAC_integrated, idents = idents_to_use,
        pt.size = 0.1, features = "Gpr133", split.by = "condition")
#not much

VlnPlot(TAC_integrated, idents = idents_to_use,
        pt.size = 0.1, features = "Emr1", cols=color_timepoint,
        split.by = "condition")
##wk5 seems a bit upregulated in MP

VlnPlot(TAC_integrated, idents = idents_to_use,
        pt.size = 0.1, features = "Cd97", cols=color_timepoint,
        split.by = "condition")
##Wk5 seems upregulated in GN, but population is low 
##a little bit same trend in MP

##Emr4 not much

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Gpr110","Gpr111", "Gpr113",
                         "Gpr115", "Gpr116"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
##everything kinda low but Gpr116

VlnPlot(TAC_integrated, idents = idents_to_use,
        pt.size = 0, features = "Gpr116", cols=color_timepoint,
        split.by = "condition")
##interestingly in CM, Gpr116 seems to go up ~5w, but drops after that
##possibly an upregulation in EC @5w too


FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Gpr56","Gpr64", "Gpr97",
                         "Gpr112", "Gpr114","Gpr126","Gpr128"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

VlnPlot(TAC_integrated, idents = idents_to_use,
        pt.size = 0.1, features = "Gpr56", cols=color_timepoint,
        split.by = "condition")
#not enough % to show violin shape, but from the dots it might be upregulated in CM at 5w

FeaturePlot(TAC_integrated, 
            reduction = "umap", 
            features = c("Lphn1","Lphn2", "Lphn3",
                         "Eltd1", "Gpr98"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

VlnPlot(TAC_integrated, idents = idents_to_use,
        pt.size = 0, features = "Eltd1", cols=color_timepoint,
        split.by = "condition")
#not that clear, but maybe higher at 8w in EC

VlnPlot(TAC_integrated, idents = idents_to_use,
        pt.size = 0.1, features = "Lphn2", cols=color_timepoint,
        split.by = "condition")
#not much


##subset##

load("data/TAC_CM_NMCC_metadata.rdata")


TAC_CM <- subset(TAC_CM_NMCC, idents = "CM") %>% 
  NormalizeData() %>% 
  ScaleData() 



##metadata$condition <- factor(metadata$condition, levels = c("0w", "2w", "5w", "8w", "11w"))

#TAC_CM@meta.data <- metadata

#once levels are included, the order of the groups are right, 
#BUT somehow all conditions will show up in heatmap, 
#even if it's filtered out and not in "cells" argument


save(TAC_CM, file = "data/TAC_CM_metadata.rdata")
###without levels 

load("data/TAC_CM_metadata.rdata")

metadata <- TAC_CM@meta.data

##set new idents (current subset only has one ident, CM)
Idents(TAC_CM) <- metadata$condition

#set levels for object; this is referring to the idents 
#(so if you don't set idents, you can't do levels right)
levels(TAC_CM) <- c("0w", "2w", "5w", "8w", "11w")


#cell_0_5w_11w <- metadata %>% 
 # rownames_to_column("cell_id") %>% 
  #filter(condition %in% c("0w", "5w", "11w"))%>% 
  #pull("cell_id")

#so originally I planned to define the cells I'm including in the heatmap, 
#but couldn't fix the order of groups 
#so gonna try subsetting out the condition that I want, and then plot

TAC_CM_plot <- subset(TAC_CM, idents = c("0w", "5w", "11w"))

#color_group <- met.brewer("Homer1", 3, direction = -1, type = "discrete")

DoHeatmap(TAC_CM_plot, 
          features = c("Myh6","Myh7", "Adrb1",
                      "Gpr116", "Gpr56", "Eltd1", "Lphn2", "Cd97", "Gpr124","Gpr133", "Emr1"),
          #group.by = "condition",
          #group.colors = color_group, #not changed in legend, a bug that's being fixed
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))


DotPlot(TAC_CM_plot,
        features = c("Myh6","Myh7", "Adrb1",
                     "Gpr116", "Gpr56", "Eltd1", "Lphn2", "Cd97", "Gpr124","Gpr133", "Emr1"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)


#####taking a brief look at FB#####
load("data/TAC_CM_NMCC_metadata_no levels.rdata")

TAC_FB <- subset(TAC_CM_NMCC, idents = "FB") %>% 
  NormalizeData() %>% 
  ScaleData() 

metadata_fb <- TAC_FB@meta.data

Idents(TAC_FB) <- metadata_fb$condition

levels(TAC_FB) <- c("0w", "2w", "5w", "8w", "11w")

DoHeatmap(TAC_FB, 
          features = c("Col3a1",
                       "Gpr116", "Gpr124","Gpr133"),
          lines.width = 20,
          #group.by = "condition",
          #group.colors = color_group, #not changed in legend, a bug that's being fixed
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))

DotPlot(TAC_FB,
        features = c("Col3a1",
                     "Gpr116", "Gpr124","Gpr133"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)
