library(tidyverse)
library(Seurat)


if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scater")

library(SeuratDisk)
library(scater)
library(reticulate)
library(hdf5r)
library(Matrix)

#downloaded data includes cell type info, BUT format in h5ad from Scanpy
#conversion from h5ad to Seurat object (not the count. including metadata and all)


Convert("raw/snRNA-seq-submission.h5ad", dest = "h5seurat", overwrite = TRUE)

obj_HDF5<- Connect("raw/snRNA-seq-submission.h5seurat", mode = "r+")

###dim discrepancy may be the issue###
x <- obj_HDF5[["assays/RNA/data"]]
sp <- sparseMatrix(i=x[["indices"]][]+1,p=x[["indptr"]][],x=x[["data"]][]  )
dim(sp)
h5attr(x = x, which = "dims")
####
obj_HDF5


human_MI <- LoadH5Seurat("raw/snRNA-seq-submission.h5seurat")
