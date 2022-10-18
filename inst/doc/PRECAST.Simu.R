## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  githubURL <- "https://github.com/feiyoung/PRECAST/blob/main/vignettes_data/data_simu.rda?raw=true"
#  download.file(githubURL,"data_simu.rda",mode='wb')
#  

## ----eval = FALSE-------------------------------------------------------------
#  load("data_simu.rda")

## ----eval = FALSE-------------------------------------------------------------
#  library(PRECAST)
#  library(Seurat)

## ----eval = FALSE-------------------------------------------------------------
#  
#  data_simu ## a list including three Seurat object with default assay: RNA

## ----eval= FALSE--------------------------------------------------------------
#  head(data_simu[[1]])

## ----eval = FALSE-------------------------------------------------------------
#  
#  ## Create
#  set.seed(2022)
#  PRECASTObj <-  CreatePRECASTObject(data_simu, customGenelist=row.names(data_simu[[1]]))
#  
#  

## ----eval = FALSE-------------------------------------------------------------
#  ## check the number of genes/features after filtering step
#  PRECASTObj@seulist
#  
#  ## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
#  PRECASTObj <-  AddAdjList(PRECASTObj, platform = "Visium")
#  
#  ## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the information in the algorithm.
#  PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal=FALSE, maxIter=30, verbose=TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  ### Given K
#  set.seed(2022)
#  PRECASTObj <- PRECAST(PRECASTObj, K=7)
#  

## ----eval = FALSE-------------------------------------------------------------
#  ## backup the fitting results in resList
#  resList <- PRECASTObj@resList
#  # PRECASTObj@resList <- resList
#  PRECASTObj <- selectModel(PRECASTObj)
#  true_cluster <- lapply(data_simu, function(x) x$true_cluster)
#  str(true_cluster)
#  mclust::adjustedRandIndex(unlist(PRECASTObj@resList$cluster), unlist(true_cluster))

## ----eval = FALSE-------------------------------------------------------------
#  
#  seuInt <- IntegrateSpaData(PRECASTObj, species='unknown')
#  seuInt
#  ## The low-dimensional embeddings obtained by PRECAST are saved in PRECAST reduction slot.

## ----eval = FALSE-------------------------------------------------------------
#  p12 <- SpaPlot(seuInt, batch=NULL,point_size=2, combine=TRUE)
#  p12
#  # users can plot each sample by setting combine=FALSE

## ----eval = FALSE-------------------------------------------------------------
#  seuInt <- AddUMAP(seuInt)
#  SpaPlot(seuInt, batch=NULL,item='RGB_UMAP',point_size=2, combine=TRUE, text_size=15)
#  
#  #seuInt <- AddTSNE(seuInt)
#  #SpaPlot(seuInt, batch=NULL,item='RGB_TSNE',point_size=2, combine=T, text_size=15)

## ----eval = FALSE-------------------------------------------------------------
#  seuInt <- AddTSNE(seuInt, n_comp = 2)
#  library(patchwork)
#  cols_cluster <- c("#E04D50", "#4374A5", "#F08A21","#2AB673", "#FCDDDE",  "#70B5B0", "#DFE0EE" ,"#D0B14C")
#  p1 <- dimPlot(seuInt,  font_family='serif', cols=cols_cluster) # Times New Roman
#  p2 <- dimPlot(seuInt, item='batch', point_size = 1,  font_family='serif')
#  p1 + p2
#  # It is noted that only sample batch 1 has cluster 4, and only sample batch 2 has cluster 7.

## ----eval = FALSE-------------------------------------------------------------
#  dimPlot(seuInt, reduction = 'UMAP3', item='cluster', cols=cols_cluster, font_family='serif')

## ----eval = FALSE-------------------------------------------------------------
#  DimPlot(seuInt, reduction = 'position')
#  DimPlot(seuInt, reduction = 'tSNE')

## ----eval = FALSE-------------------------------------------------------------
#  dat_deg <- FindAllMarkers(seuInt)
#  library(dplyr)
#  n <- 2
#  dat_deg %>%
#    group_by(cluster) %>%
#    top_n(n = n, wt = avg_log2FC) -> top10
#  
#  head(top10)
#  

## ----eval = FALSE-------------------------------------------------------------
#  sessionInfo()

