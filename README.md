# scRNA-seq
 
This is a series projects involving published scRNA-seq datasets for hearts in mice and humans in normal vs diseased states. The goal was to re-run the analysis and visualization pipeline using Seurat focusing on the genes of interest in our lab, as well as to build Shiny dashboards for the lab to explore the data. 

## data/
In this folder, each folder contains a dataset that I re-run analysis for (raw data, clean data objects, R script). Here are the cool papers that published the datasets that I used: 
 - [Tucker, Nathan R., et al. "Transcriptional and cellular diversity of the human heart."Circulation 142.5 (2020): 466-482.](https://doi.org/10.1161/CIRCULATIONAHA.119.045401)
 - [Kuppe, C., Ramirez Flores, R.O., Li, Z. et al. Spatial multi-omic map of human myocardial infarction. Nature 608, 766–777 (2022).](https://doi.org/10.1038/s41586-022-05060-x)
 - [Wang, L., Yu, P., Zhou, B. et al. Single-cell reconstruction of the adult human heart during heart failure and recovery reveals the cellular landscape underlying cardiac function. Nat Cell Biol 22, 108–119 (2020).](https://doi.org/10.1038/s41556-019-0446-7)
 - [Litviňuková, M., Talavera-López, C., Maatz, H. et al. Cells of the adult human heart. Nature 588, 466–472 (2020).](https://doi.org/10.1038/s41586-020-2797-4)
 - [Rao, M., Wang, X., Guo, G. et al. Resolving the intertwining of inflammation and fibrosis in human heart failure at single-cell level. Basic Res Cardiol 116, 55 (2021).](https://doi.org/10.1007/s00395-021-00897-1)
 - [Koenig, A.L., Shchukina, I., Amrute, J. et al. Single-cell transcriptomics reveals cell-type-specific diversification in human heart failure. Nat Cardiovasc Res 1, 263–280 (2022).](https://doi.org/10.1038/s44161-022-00028-6)
 - [Ren, Zongna, et al. "Single-cell reconstruction of progression trajectory reveals intervention principles in pathological cardiac hypertrophy." Circulation 141.21 (2020): 1704-1719.](https://doi.org/10.1161/CIRCULATIONAHA.119.043053)
 - [Farbehi, Nona, et al. "Single-cell expression profiling reveals dynamic flux of cardiac stromal, vascular and immune cells in health and injury." Elife 8 (2019): e43882.](https://doi.org/10.7554/eLife.43882)
 - [Martini, Elisa, et al. "Single-cell sequencing of mouse heart immune infiltrate in pressure overload–driven heart failure reveals extent of immune activation." Circulation 140.25 (2019): 2089-2107.](https://doi.org/10.1161/CIRCULATIONAHA.119.041694)
 - [Revelo, Xavier S., et al. "Cardiac resident macrophages prevent fibrosis and stimulate angiogenesis." Circulation research 129.12 (2021): 1086-1101.](https://doi.org/10.1161/CIRCRESAHA.121.319737)

The data R project deals with clean data objects from multiple datasets and preps data for shiny apps later on.
The shiny app_multiple dataset backup is a backup for record on how to update SelectInput based on another input (observeEvent). This was not used in shiny apps because the datasets were too big to put into one app, and the update input feature was not needed for single dataset.

## quarto/
This is to put together a report looking at specific genes (adhesion GPCRs) in lymphatic endothelial cells in two mouse heart scRNA-seq datasets. I used quarto to try it out, but not really doing justice to the new amazing functionality of quarto as I'm kinda just using it as Rmarkdown with html output. But fun..!

## shiny-xx/
The three shiny folders are the shiny apps for the lab to explore the gene expression difference across different conditions in major cell types of the heart (with dotplot and violin plot), with one dataset for each app. Similar design, with small modifications on cell type choices in some datasets. 



