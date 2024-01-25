## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy = TRUE,
  fig.width = 11,
  tidy.opts = list(width.cutoff = 95),
  message = FALSE,
  warning = FALSE,
  time_it = TRUE
)
all_times <- list()  # store the time for each chunk
knitr::knit_hooks$set(time_it = local({
  now <- NULL
  function(before, options) {
    if (before) {
      now <<- Sys.time()
    } else {
      res <- difftime(Sys.time(), now, units = "secs")
      all_times[[options$label]] <<- res
    }
  }
}))

## ----eval =FALSE--------------------------------------------------------------
#  suppressPackageStartupMessages(library(Seurat))
#  suppressPackageStartupMessages(library(SingleCellExperiment))
#  name_ID4 <- as.character(c(151673, 151674, 151675, 151676))
#  
#  ### Read data in an online manner
#  n_ID <- length(name_ID4)
#  url_brainA <- "https://github.com/feiyoung/DR-SC.Analysis/raw/main/data/DLPFC_data/"; url_brainB <- ".rds"
#  seuList <- list()
#  if(!require(ProFAST)){
#    remotes::install_github("feiyoung/ProFAST")
#  }
#  for(i in 1:n_ID){
#    # i <- 1
#    cat('input brain data', i, '\n')
#    # load and read data
#    dlpfc <- readRDS(url(paste0(url_brainA, name_ID4[i],url_brainB) ))
#    count <- dlpfc@assays@data$counts
#    row.names(count) <- ProFAST::transferGeneNames(row.names(count), species = "Human")
#    seu1 <- CreateSeuratObject(counts = count,
#                               meta.data = as.data.frame(colData(dlpfc)),
#                                min.cells = 10,min.features = 10)
#    seuList[[i]] <- seu1
#  }
#  # saveRDS(seuList, file='seuList4.RDS')
#  

## ----eval = FALSE-------------------------------------------------------------
#  library(PRECAST)

## ----eval = FALSE-------------------------------------------------------------
#  seuList <- readRDS("seuList4.RDS")
#  seuList ## a list including  Seurat objects

## ----eval = FALSE-------------------------------------------------------------
#  metadataList <- lapply(seuList, function(x) x@meta.data)
#  
#  for(r in seq_along(metadataList)){
#    meta_data <- metadataList[[r]]
#    cat(all(c("row", "col") %in% colnames(meta_data)), '\n') ## the names are correct!
#  
#  }

## ----eval = FALSE-------------------------------------------------------------
#  set.seed(2023)
#  preobj <- CreatePRECASTObject(seuList = seuList, selectGenesMethod="HVGs",
#                                gene.number = 2000) #

## ----eval = FALSE-------------------------------------------------------------
#  ## check the number of genes/features after filtering step
#  preobj@seulist
#  ## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
#  PRECASTObj <-  AddAdjList(preobj, platform = "Visium")
#  ## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the
#  ## information in the algorithm.
#  PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = TRUE, coreNum = 1,
#                              maxIter=30, verbose = TRUE)
#  

## ----eval = FALSE-------------------------------------------------------------
#  ### Given K
#  PRECASTObj <- PRECAST(PRECASTObj, K= 7)

## ----eval = FALSE-------------------------------------------------------------
#  ## backup the fitting results in resList
#  resList <- PRECASTObj@resList
#  PRECASTObj <- SelectModel(PRECASTObj)
#  ari_precast <- sapply(1:length(seuList), function(r) mclust::adjustedRandIndex(PRECASTObj@resList$cluster[[r]], PRECASTObj@seulist[[r]]$layer_guess_reordered))
#  mat <- matrix(round(ari_precast,2), nrow=1)
#  name_ID4 <- as.character(c(151673, 151674, 151675, 151676))
#  colnames(mat) <-   name_ID4
#  DT::datatable(mat)

## ----eval =FALSE--------------------------------------------------------------
#  PRECASTObj2 <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 4,
#                              maxIter=30, verbose = TRUE) # set 4 cores to run in parallel.
#  PRECASTObj2 <- PRECAST(PRECASTObj2, K= 5:8)
#  ## backup the fitting results in resList
#  resList2 <- PRECASTObj2@resList
#  PRECASTObj2 <- SelectModel(PRECASTObj2)
#  

## ----eval = FALSE-------------------------------------------------------------
#  print(PRECASTObj@seuList)
#  seuInt <- IntegrateSpaData(PRECASTObj, species='Human')
#  seuInt
#  ## The low-dimensional embeddings obtained by PRECAST are saved in PRECAST reduction slot.

## ----eval =FALSE--------------------------------------------------------------
#  ## assign the raw Seurat list object to it.
#  ## For illustration, we generate a new seuList with more genes;
#  ## For integrating all genes, users can set `seuList <- bc2`.
#  genes <- c(row.names(PRECASTObj@seulist[[1]]), row.names(seuList[[1]])[1:10])
#  seuList_sub <- lapply(seuList, function(x) x[genes,])
#  PRECASTObj@seuList <- seuList_sub #
#  seuInt <- IntegrateSpaData(PRECASTObj, species='Human')
#  seuInt

## ----eval =FALSE--------------------------------------------------------------
#  PRECASTObj@seuList <- NULL
#  ## At the same time, we can set subsampling to speed up the computation.
#  seuInt <- IntegrateSpaData(PRECASTObj, species='Human', seuList=seuList_sub, subsample_rate = 0.5)
#  seuInt

## ----eval = FALSE-------------------------------------------------------------
#  cols_cluster <- chooseColors(palettes_name = 'Classic 20', n_colors=7, plot_colors = TRUE)

## ----eval = FALSE, fig.height = 8, fig.width=9--------------------------------
#  
#  p12 <- SpaPlot(seuInt, item='cluster', batch=NULL,point_size=1, cols=cols_cluster, combine=TRUE, nrow.legend=7)
#  p12
#  
#  # users can plot each sample by setting combine=FALSE

## ----eval = FALSE, fig.height = 8, fig.width=8.5------------------------------
#  library(ggplot2)
#  pList <- SpaPlot(seuInt, item='cluster', batch=NULL,point_size=2.5, cols=cols_cluster, combine=FALSE, nrow.legend=7)
#  pList <- lapply(pList, function(x) x +  coord_flip() + scale_x_reverse())
#  drawFigs(pList, layout.dim = c(2,2), common.legend = TRUE, legend.position = 'right', align='hv')
#  

## ----eval = FALSE, fig.height = 6, fig.width=5.5------------------------------
#  seuInt <- AddUMAP(seuInt)
#  p13List <- SpaPlot(seuInt, batch=NULL,item='RGB_UMAP',point_size=2, combine=FALSE, text_size=15)
#  p13List <- lapply(p13List, function(x) x +  coord_flip() + scale_x_reverse())
#  drawFigs(p13List, layout.dim = c(2,2), common.legend = TRUE, legend.position = 'right', align='hv')
#  #seuInt <- AddTSNE(seuInt)
#  #SpaPlot(seuInt, batch=NULL,item='RGB_TSNE',point_size=2, combine=T, text_size=15)

## ----eval = FALSE, fig.height = 4.5, fig.width=12-----------------------------
#  seuInt <- AddTSNE(seuInt, n_comp = 2)
#  p1 <- dimPlot(seuInt, item='cluster', point_size = 0.5, font_family='serif', cols=cols_cluster,border_col="gray10", nrow.legend=14, legend_pos='right') # Times New Roman
#  p2 <- dimPlot(seuInt, item='batch', point_size = 0.5,  font_family='serif', legend_pos='right')
#  
#  drawFigs(list(p1, p2), layout.dim = c(1,2), legend.position = 'right', align='hv')
#  

## ----eval = FALSE-------------------------------------------------------------
#  library(Seurat)
#  dat_deg <- FindAllMarkers(seuInt)
#  library(dplyr)
#  n <- 5
#  dat_deg %>%
#    group_by(cluster) %>%
#    top_n(n = n, wt = avg_log2FC) -> top10
#  
#  

## ----eval = FALSE, fig.height = 6, fig.width=8--------------------------------
#  library(ggplot2)
#  ## HeatMap
#  p1 <- DotPlot(seuInt, features = unique(top10$gene), col.min = 0, col.max = 1) +
#    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8))
#  p1

## -----------------------------------------------------------------------------
sessionInfo()

