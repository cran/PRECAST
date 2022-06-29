## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  library(Seurat)
#  library(PRECAST)
#  seuList <- gendata_seulist(height1=20, width1=20, height2=24, width2=25,p=200, K=4, alpha=20, sigma2=0.1)
#  seuList
#  head(seuList[[1]])
#  ## Must include the columns named "row" and "col" for saving the spatial coordinates
#  

## ----eval = FALSE-------------------------------------------------------------
#  ## Create
#  PRECASTObj <-  CreatePRECASTObject(seuList)
#  

## ----eval = FALSE-------------------------------------------------------------
#  ## check the number of genes/features after filtering step
#  PRECASTObj@seulist
#  
#  ## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
#  PRECASTObj <-  AddAdjList(PRECASTObj, platform = "ST")
#  
#  ## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the information in the algorithm.
#  PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal=TRUE, verbose=TRUE, seed=2022)

## ----eval = FALSE-------------------------------------------------------------
#  ### Given K
#  PRECASTObj <- PRECAST(PRECASTObj, K=5)
#  

## ----eval = FALSE-------------------------------------------------------------
#  ## backup the fitting results in resList
#  resList <- PRECASTObj@resList
#  # PRECASTObj@resList <- resList
#  PRECASTObj <- selectModel(PRECASTObj)
#  true_cluster <- lapply(seuList, function(x) x$true_cluster)
#  str(true_cluster)
#  mclust::adjustedRandIndex(unlist(PRECASTObj@resList$cluster), unlist(true_cluster))

## ----eval = FALSE-------------------------------------------------------------
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
#  DimPlot(seuInt, reduction = 'PRECAST')

## ----eval = FALSE-------------------------------------------------------------
#  dat_deg <- FindAllMarkers(seuInt)
#  library(dplyr)
#  n <- 10
#  dat_deg %>%
#    group_by(cluster) %>%
#    top_n(n = n, wt = avg_log2FC) -> top10
#  
#  seuInt <- ScaleData(seuInt)
#  seus <- subset(seuInt, downsample = 400)
#  color_id <- as.numeric(levels(Idents(seus)))
#  
#  library(ggplot2)
#  
#  ## HeatMap
#  p1 <- doHeatmap(seus, features = top10$gene, cell_label= "Domain",
#                  grp_label = F, grp_color = cols_cluster,
#                  pt_size=6,slot = 'scale.data') +
#    theme(legend.text = element_text(size=16),
#          legend.title = element_text(size=18, face='bold'),
#          axis.text.y = element_text(size=7, face= "italic", family='serif'))
#  p1
#  

## ----eval = FALSE-------------------------------------------------------------
#  sessionInfo()
