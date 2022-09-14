#' Data Prep

library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)
library(data.table)


###import
######Fread() handles big files!!!!!
counts <- fread("raw/GSE183852_Integrated_Counts.csv")
##took 2min or longer, 23g file
#human_DCM <- CreateSeuratObject(counts, project = "human DCM")
##to big, can't be processed


## 1. Test the split and create seurate then merge
## Splip and process
df1 = counts[,c(1,10:13)]


## 2. Partition Table into five parts

### Split counts into seperate files
dim(counts)
## there are 45,068 rows. 269795 columns 
#cannot partition by rows, because later on merging objects requires unique cell id
#so have to do it by columns, while keeping first column (gene names)
counts[,1:35000] %>% fwrite('processed data/GSE183852_Integrated_Counts_pt1.csv')
counts[,c(1, 35001:70000)] %>% fwrite('processed data/GSE183852_Integrated_Counts_pt2.csv')
counts[,c(1, 70001:105000)] %>% fwrite('processed data/GSE183852_Integrated_Counts_pt3.csv')
counts[, c(1, 105001:140000)] %>% fwrite('processed data/GSE183852_Integrated_Counts_pt4.csv')
counts[,c(1, 140001:175000)] %>% fwrite('processed data/GSE183852_Integrated_Counts_pt5.csv')
counts[,c(1, 175001:210000)] %>% fwrite('processed data/GSE183852_Integrated_Counts_pt6.csv')
counts[, c(1, 210001:245000)] %>% fwrite('processed data/GSE183852_Integrated_Counts_pt7.csv')
counts[, c(1, 245001:269795)] %>% fwrite('processed data/GSE183852_Integrated_Counts_pt8.csv')
