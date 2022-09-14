# 0. Setup -----
library(tidyverse)

# 1. Prepare UI object -----
load("healthy MI non myocyte cardiac/clean/all_features_MI.rdata")

features_MI <- tibble(gene = all_features, dataset = "2019 eLife mouse MI")
 
load("GSE120064_TAC CM_NMCC_mtx_cell_info/clean/all_features_TAC.rdata")  
  
features_TAC <- tibble(gene = all_features_TAC, dataset = "2020 Circulation mouse TAC")

ui_features <- bind_rows(features_MI, features_TAC)
  
save(ui_features, file = "../shiny/R/ui_features.rdata")

# 2. Prepare Server object ----

load("../data/healthy MI non myocyte cardiac/clean/MI_EC_normalized.rdata") #object name: TIP_MI_EC
load("../data/healthy MI non myocyte cardiac/clean/MI_MP_normalized.rdata") #TIP_MI_MP
load("../data/GSE120064_TAC CM_NMCC_mtx_cell_info/clean/TAC_EC_metadata.rdata") #TAC_EC
load("../data/GSE120064_TAC CM_NMCC_mtx_cell_info/clean/TAC_MP_metadata.rdata") #TAC_MP

server_data = tibble(
  dataset = c("2019 eLife mouse MI", "2020 Circulation mouse TAC"),
  EC = list(TIP_MI_EC, TAC_EC),
  MP = list(TIP_MI_MP,TAC_MP)
) %>% 
  pivot_longer(c(EC,MP), names_to = 'cell', values_to = 'obj')


save(server_data, file = "../shiny/R/server_data.rdata")
