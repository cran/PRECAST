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
#  data_simu ## a list including three Seurat object with default assay: RNA

## ----eval= FALSE--------------------------------------------------------------
#  head(data_simu[[1]])

## ----eval= FALSE--------------------------------------------------------------
#  row.names(data_simu[[1]])[1:10]

## ----eval = FALSE-------------------------------------------------------------
#  ## Get the gene-by-spot read count matrices
#  countList <- lapply(data_simu, function(x){
#    assay <- DefaultAssay(x)
#    GetAssayData(x, assay = assay, slot='counts')
#  
#  } )
#  
#  ## Check the spatial coordinates: Yes, they are named as "row" and "col"!
#  head(data_simu[[1]]@meta.data)
#  
#  ## Get the meta data of each spot for each data batch
#  metadataList <- lapply(data_simu, function(x) x@meta.data)
#  
#  
#  ## ensure the row.names of metadata in metaList are the same as that of colnames count matrix in countList
#  M <- length(countList)
#  for(r in 1:M){
#    row.names(metadataList[[r]]) <- colnames(countList[[r]])
#  }
#  
#  
#  ## Create the Seurat list  object
#  
#  seuList <- list()
#  for(r in 1:M){
#    seuList[[r]] <- CreateSeuratObject(counts = countList[[r]], meta.data=metadataList[[r]], project = "PRECASTsimu")
#  }
#  

## ----eval = FALSE-------------------------------------------------------------
#  
#  ## Create PRECASTObject
#  set.seed(2022)
#  PRECASTObj <-  CreatePRECASTObject(seuList, customGenelist=row.names(seuList[[1]]))
#  
#  ## User can retain the raw seuList by the following commond.
#  ##  PRECASTObj <-  CreatePRECASTObject(seuList, customGenelist=row.names(seuList[[1]]), rawData.preserve = TRUE)
#  

## ----eval = FALSE-------------------------------------------------------------
#  ## check the number of genes/features after filtering step
#  PRECASTObj@seulist
#  
#  ## seuList is null since the default value `rawData.preserve` is FALSE.
#  PRECASTObj@seuList
#  
#  ## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
#  PRECASTObj <-  AddAdjList(PRECASTObj, platform = "Visium")
#  
#  ## Add a model setting in advance for a PRECASTObj object: verbose =TRUE helps outputing the information in the algorithm; coreNum set the how many cores are used in PRECAST. If you run PRECAST for multiple number of clusters, you can set multiple cores; otherwise, set it to 1.
#  PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal=FALSE, maxIter=30, verbose=TRUE,
#                               coreNum =1)

## ----eval = FALSE-------------------------------------------------------------
#  ### Given K
#  set.seed(2022)
#  PRECASTObj <- PRECAST(PRECASTObj, K=7)

## ----eval =FALSE--------------------------------------------------------------
#  ## Reset  parameters by increasing cores.
#  PRECASTObj2 <- AddParSetting(PRECASTObj, Sigma_equal=FALSE, maxIter=30, verbose=TRUE,
#                               coreNum =2)
#  set.seed(2023)
#  PRECASTObj2 <- PRECAST(PRECASTObj2, K=6:7)
#  
#  resList2 <- PRECASTObj2@resList
#  PRECASTObj2 <- SelectModel(PRECASTObj2)
#  

## ----eval = FALSE-------------------------------------------------------------
#  ## check the fitted results: there are four list for the fitted results of each K (6:9).
#  str(PRECASTObj@resList)
#  ## backup the fitted results in resList
#  resList <- PRECASTObj@resList
#  # PRECASTObj@resList <- resList
#  PRECASTObj <- SelectModel(PRECASTObj)
#  ## check the best and re-organized results
#  str(PRECASTObj@resList) ## The selected best K is 7

## ----eval = FALSE-------------------------------------------------------------
#  true_cluster <- lapply(PRECASTObj@seulist, function(x) x$true_cluster)
#  str(true_cluster)
#  mclust::adjustedRandIndex(unlist(PRECASTObj@resList$cluster), unlist(true_cluster))

## ----eval = FALSE-------------------------------------------------------------
#  
#  seuInt <- IntegrateSpaData(PRECASTObj, species='unknown')
#  seuInt
#  ## The low-dimensional embeddings obtained by PRECAST are saved in PRECAST reduction slot.

## ----eval = FALSE-------------------------------------------------------------
#  cols_cluster <- chooseColors(palettes_name = 'Nature 10', n_colors = 7, plot_colors = TRUE)

## ----eval = FALSE, fig.height=5, fig.width=7----------------------------------
#  p12 <- SpaPlot(seuInt, batch=NULL, cols=cols_cluster, point_size=2, combine=TRUE)
#  p12
#  # users can plot each sample by setting combine=FALSE

## ----eval = FALSE, fig.height=2.6, fig.width=7--------------------------------
#  pList <- SpaPlot(seuInt, batch=NULL, cols=cols_cluster, point_size=2, combine=FALSE, title_name=NULL)
#  drawFigs(pList[1:2], layout.dim = c(1,2), common.legend = TRUE, legend.position = 'right', align='hv')
#  

## ----eval = FALSE, fig.height=5, fig.width=6----------------------------------
#  seuInt <- AddUMAP(seuInt)
#  SpaPlot(seuInt, batch=NULL,item='RGB_UMAP',point_size=1, combine=TRUE, text_size=15)
#  
#  ## Plot tSNE RGB plot
#  #seuInt <- AddTSNE(seuInt)
#  #SpaPlot(seuInt, batch=NULL,item='RGB_TSNE',point_size=2, combine=T, text_size=15)

## ----eval = FALSE, fig.height=8, fig.width=6----------------------------------
#  seuInt <- AddTSNE(seuInt, n_comp = 2)
#  
#  p1 <- dimPlot(seuInt, item='cluster', font_family='serif', cols=cols_cluster) # Times New Roman
#  p2 <- dimPlot(seuInt, item='batch', point_size = 1,  font_family='serif')
#  drawFigs(list(p1, p2), common.legend=FALSE, align='hv')
#  # It is noted that only sample batch 1 has cluster 4, and only sample batch 2 has cluster 7.

## ----eval = FALSE, fig.height=4, fig.width=6----------------------------------
#  dimPlot(seuInt, reduction = 'UMAP3', item='cluster', cols=cols_cluster, font_family='serif')

## ----eval = FALSE, fig.height=3, fig.width=8----------------------------------
#  library(Seurat)
#  p1 <- DimPlot(seuInt[,1: 4226], reduction = 'position', cols=cols_cluster, pt.size =1) # plot the first data batch: first 4226 spots.
#  p2 <- DimPlot(seuInt, reduction = 'tSNE',cols=cols_cluster, pt.size=1)
#  drawFigs(list(p1, p2), layout.dim = c(1,2), common.legend = TRUE)

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

## -----------------------------------------------------------------------------
sessionInfo()

