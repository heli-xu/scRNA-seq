# 0. Setup -----
library(tidyverse)

# 1. Prepare UI object -----
{
  load("healthy MI non myocyte cardiac/clean/all_features_MI.rdata")
  
  features_MI <- tibble(gene = all_features, dataset = "2019 eLife mouse MI")
  
  load("GSE120064_TAC CM_NMCC_mtx_cell_info/clean/all_features_TAC.rdata")  
  
  features_TAC <- tibble(gene = all_features_TAC, dataset = "2020 Circulation mouse TAC")
  
  ui_features_all <- bind_rows(features_MI, features_TAC)
}
##for individual dataset, all_features files were copied directly

# 2. Prepare Server object ----(2 mice + 1 human datasets)

load("../data/healthy MI non myocyte cardiac/clean/MI_EC_diet.rdata") #object name same as file name
load("../data/healthy MI non myocyte cardiac/clean/MI_MP_diet.rdata") 
load("../data/GSE120064_TAC CM_NMCC_mtx_cell_info/clean/TAC_EC_diet.rdata") 
load("../data/GSE120064_TAC CM_NMCC_mtx_cell_info/clean/TAC_MP_diet.rdata") 
load("../data/GSE120064_TAC CM_NMCC_mtx_cell_info/clean/TAC_CM_diet.rdata")
load("../data/2020 nat cell biol human normal vs failing heart/clean/human_CM_diet.rdata")
load("../data/2020 nat cell biol human normal vs failing heart/clean/human_EC_diet.rdata")
load("../data/2020 nat cell biol human normal vs failing heart/clean/human_MP_diet.rdata")

server_data_all = tibble(
  dataset = c("2019 eLife mouse MI", "2020 Circulation mouse TAC", "2020 NatCellBio human HF"),
  EC = list(MI_EC_diet, TAC_EC_diet, human_EC_diet),
  MP = list(MI_MP_diet,TAC_MP_diet, human_MP_diet),
  CM = list(NA, TAC_CM_diet, human_CM_diet)
) %>% 
  pivot_longer(c(EC,MP,CM), names_to = 'cell', values_to = 'obj')


pull_dataset <- function(x) {
  filter(server_data_all, dataset == x) %>% 
    select(-dataset)
}

server_data_MI <- pull_dataset("2019 eLife mouse MI")
save(server_data_MI, file = "../shiny-mouseMI/R/server_data_MI.rdata")

server_data_TAC <- pull_dataset("2020 Circulation mouse TAC")

save(server_data_TAC, file = "../shiny-mouseTAC/R/server_data_TAC.rdata")

server_data_humanHF <- pull_dataset("2020 NatCellBio human HF")
save(server_data_humanHF, file = "../shiny-humanHF/R/server_data_humanHF.rdata")
