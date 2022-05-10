library(dplyr)
library(Seurat)
library(Matrix)
library(tidyverse)

#########import

pool1.data <- Read10X(data.dir = "raw/pool1")
pool2.data <- Read10X(data.dir = "raw/pool2")
