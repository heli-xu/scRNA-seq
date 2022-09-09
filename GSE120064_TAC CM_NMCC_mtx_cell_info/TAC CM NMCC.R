library(dplyr)
library(Seurat)

library(tidyverse)
library(data.table)
library(MetBrewer)
library(ggplot2)
library(viridis)

###load raw data and cell info (metadata)##

counts <- fread("raw/TAC_raw_umi_matrix.csv")

rownames(counts) <- counts$V1
#for matrix, rownames_to_column doesnt seem to work well

counts$V1 <- NULL
##800m, doable, but stuck at next step

TAC_CM_NMCC <- CreateSeuratObject(counts, project = "TAC CM NMCC",
                                  min.cells =3, min.features  = 200)
save(TAC_CM_NMCC, file = "data/TAC_CM_NMCC.rdata")

##too big, although the object turns out to be 400m

load("data/TAC_CM_NMCC.rdata")

cell_info <- read.table("raw/GSE120064_TAC_clean_cell_info_summary.txt", header = T)

####metadata wrangling####

metadata <- TAC_CM_NMCC@meta.data %>% 
  rownames_to_column(var = "cell_id")

cell_info <- cell_info %>% 
  rename(cell_id = CellID)

metadata_full <- metadata %>% 
  left_join(cell_info, by = "cell_id") %>% 
  column_to_rownames("cell_id")
##keep the rownames as cell_id

##metadata_full$condition <- factor(metadata_full$condition, levels = c("0w", "2w","5w","8w","11w"))
##good extra thing to do for future visualization

TAC_CM_NMCC@meta.data <- metadata_full

Idents(TAC_CM_NMCC) <- metadata_full$CellType

save(TAC_CM_NMCC, file = "data/TAC_CM_NMCC_metadata_celltype.rdata")
##conditions no levels, cell type added


####Normalization without integration
TAC_no_integ <- NormalizeData(TAC_CM_NMCC) %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters(resolution = 1) %>% 
  RunUMAP(dim=1:20)

DimPlot(TAC_no_integ, reduction = "umap", group.by = "condition")
##some cluster seems only to show one condition, but paper used some housekeeping genes
#to demonstrate minimal batch effects


####Sctransform with integration ###

load("data/cell_cycle_mouse.rdata")

s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")


split_TAC <- SplitObject(TAC_CM_NMCC, split.by = "condition")

split_TAC <- split_TAC[c("0w","2w","5w","8w","11w")] %>%  
  map(~.x %>% 
        NormalizeData() %>% 
        CellCycleScoring(g2m.features = g2m_genes, 
                         s.features = s_genes) %>% 
        SCTransform(verbose=FALSE))
##took quite a few min to run

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
##took quite long (15min for 8G RAM)


# Integrate across conditions
TAC_integrated <- IntegrateData(anchorset = integ_anchors,
                                normalization.method = "SCT",
                                verbose=FALSE)
##cannot run with 8G RAM


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
metadata %>% 
  ggplot(aes(x=condition, fill=CellType)) +
  geom_bar() +
  scale_fill_manual(values = met.brewer("Renoir",8, direction = -1, type = "discrete"))+
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


####subset CM####

load("data/TAC_CM_NMCC_metadata_celltype.rdata")


TAC_CM <- subset(TAC_CM_NMCC, idents = "CM") %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>%
  ScaleData() 

##FindVariableFeature significantly reduced TAC_CM size
#for heatmap, if gene not in variable feature, will not show up, so better to include all features;
#for dotplot and violin plot, it's fine-will be very small anyway

##metadata$condition <- factor(metadata$condition, levels = c("0w", "2w", "5w", "8w", "11w"))

#TAC_CM@meta.data <- metadata

#once levels are included, the order of the groups are right, 
#BUT somehow all conditions will show up in heatmap, 
#even if it's filtered out and not in "cells" argument



metadata_cm <- TAC_CM@meta.data

##set new idents (current subset only has one ident, CM)
Idents(TAC_CM) <- metadata_cm$condition

#set levels for object; this is referring to the idents 
#(so if you don't set idents, you can't do levels right)
levels(TAC_CM) <- c("0w", "2w", "5w", "8w", "11w")


save(TAC_CM, file = "data/TAC_CM_metadata.rdata")
###with levels as 0w, 2w, 5w, 8w, 11w
##scaledata with variable features


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

DoHeatmap(TAC_CM, 
          features = c("Myh6","Myh7", "Gpr116"),
          #group.by = "condition",
          #group.colors = color_group, #not changed in legend, a bug that's being fixed
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))

DotPlot(TAC_CM,
        features = c("Myh6","Myh7", "Adrb1",
                     "Gpr116", "Gpr56", "Eltd1", "Lphn2", "Cd97", "Gpr124","Gpr133", "Emr1"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

DotPlot(TAC_CM_plot,
        features = c("Myh7","Gpr116"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


VlnPlot(TAC_CM, features = "Gpr116")

##saving all features as rdata##
#check if row number match the mother object (TAC_CM_NMCC)

all_features_TAC <- rownames(TAC_CM@assays$RNA@data)

save(all_features_TAC, file = "data/all_features_TAC.rdata")


#####taking a brief look at FB#####
load("data/TAC_CM_NMCC_metadata_celltype.rdata")

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

#####ADGRs in EC####
load("data/TAC_CM_NMCC_metadata_celltype.rdata")
  
TAC_EC <- subset(TAC_CM_NMCC, idents = "EC") %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>%
  NormalizeData() %>% 
  ScaleData() 

metadata_ec <- TAC_EC@meta.data

Idents(TAC_EC) <- metadata_ec$condition

levels(TAC_EC) <- c("0w", "2w", "5w", "8w", "11w")

save(TAC_EC, file = "data/TAC_EC_metadata.rdata")

DoHeatmap(TAC_EC, 
          features = c("Pecam1",
                       "Gpr116", "Gpr56", "Eltd1", "Lphn2", 
                       "Cd97", "Gpr124","Gpr133", "Emr1"),
          lines.width = 20,
          #group.by = "condition",
          #group.colors = color_group, #not changed in legend, a bug that's being fixed
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))

DotPlot(TAC_EC,
        features = c("Pecam1",
                     "Gpr116", "Gpr56", "Eltd1", "Lphn2", 
                     "Cd97", "Gpr124","Gpr133", "Emr1"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

VlnPlot(TAC_EC, features = "Gpr116")
VlnPlot(TAC_EC, features = "Eltd1")

#####ADGRs in MP####
load("data/TAC_CM_NMCC_metadata_celltype.rdata")
#with celltype info

metadata <- TAC_CM_NMCC@meta.data

TAC_MP <- subset(TAC_CM_NMCC, idents = "MP") %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, verbose = FALSE) %>%
  ScaleData() 

metadata_mp <- TAC_MP@meta.data

Idents(TAC_MP) <- metadata_mp$condition

levels(TAC_MP) <- c("0w", "2w", "5w", "8w", "11w")

save(TAC_MP, file = "data/TAC_MP_metadata.rdata")

DoHeatmap(TAC_MP, 
          features = c("Cd14",
                       "Gpr116", "Gpr56",
                       "Emr1","Emr2", "Emr3",
                       "Emr4", "Cd97",
                       "Gpr124","Gpr133", "Lphn2","Eltd1"),
          lines.width = 20,
          #group.by = "condition",
          #group.colors = color_group, #not changed in legend, a bug that's being fixed
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))

DotPlot(TAC_MP,
        features = c("Cd14",
                     "Gpr116", "Gpr56",
                     "Emr1","Emr2", "Emr3",
                     "Emr4", "Cd97",
                     "Gpr124","Gpr133", "Lphn2","Eltd1"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

VlnPlot(TAC_MP, features = c("Gpr116","Emr1"))

#####ADGRs in T####
TAC_T <- subset(TAC_CM_NMCC, idents = "T") %>% 
  NormalizeData() %>% 
  ScaleData() 

metadata_t <- TAC_T@meta.data

Idents(TAC_T) <- metadata_t$condition

levels(TAC_T) <- c("0w", "2w", "5w", "8w", "11w")

DoHeatmap(TAC_T, 
          features = c("Cd3d",
                       "Gpr116", "Gpr56",
                       "Emr1","Emr2", "Emr3",
                       "Emr4", "Cd97",
                       "Gpr124","Gpr133", "Lphn2","Eltd1"),
          lines.width = 20,
          #group.by = "condition",
          #group.colors = color_group, #not changed in legend, a bug that's being fixed
          disp.min = -0.5, disp.max = 3)+
  scale_fill_viridis(option = "B", na.value = "white")+ #na.value for white line between groups
  theme(text = element_text(size = 20))

DotPlot(TAC_T,
        features = c(#"Cd3d",
                     "Gpr116", "Gpr56",
                     "Emr1","Emr2", "Emr3",
                     "Emr4", "Cd97",
                     "Gpr124","Gpr133", "Lphn2","Eltd1"),
        col.min = 0, col.max = 3,
        dot.min = 0, dot.scale = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

VlnPlot(TAC_T, features = "Gpr116")
#"Gpr56","Cd97" no violin shape
