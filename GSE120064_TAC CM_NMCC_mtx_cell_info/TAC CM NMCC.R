library(dplyr)
library(Seurat)

library(tidyverse)
library(data.table)
library(MetBrewer)


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


TAC_CM_NMCC@meta.data <- metadata_full

save(TAC_CM_NMCC, file = "data/TAC_CM_NMCC_metadata.rdata")

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
  theme_minimal() +
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

#VlnPlot(TAC_integrated, pt.size = 0, features = "Gpr124")
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

#VlnPlot(TAC_integrated, pt.size = 0, features = "Gpr133")
#not much

VlnPlot(TAC_integrated, idents = idents_to_use,
        pt.size = 0, features = "Emr1", cols=color_timepoint,
        split.by = "condition")
##wk5 seems a bit upregulated in MP

VlnPlot(TAC_integrated, idents = idents_to_use,
        pt.size = 0, features = "Cd97", cols=color_timepoint,
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