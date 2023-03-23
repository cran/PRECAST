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

## ----eval=FALSE---------------------------------------------------------------
#  githubURL <- "https://github.com/feiyoung/PRECAST/blob/main/vignettes_data/dlpfc_151672.rda?raw=true"
#  download.file(githubURL,"dlpfc_151672.rda",mode='wb')
#  

## ----eval = FALSE-------------------------------------------------------------
#  load("dlpfc_151672.rda")

## ----eval = FALSE-------------------------------------------------------------
#  library(PRECAST)
#  library(Seurat)

## ----eval = FALSE-------------------------------------------------------------
#  dlpfc_151672 ## a list including two Seurat object

## ----eval = FALSE-------------------------------------------------------------
#  library(PRECAST)
#  preobj <- CreatePRECASTObject(seuList = list(dlpfc_151672))

## ----eval = FALSE-------------------------------------------------------------
#  ## check the number of genes/features after filtering step
#  preobj@seulist
#  ## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
#  PRECASTObj <-  AddAdjList(preobj, platform = "Visium")
#  ## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the
#  ## information in the algorithm.
#  PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 1,
#                              maxIter=30, verbose = TRUE)
#  

## ----eval = FALSE-------------------------------------------------------------
#  ### Given K
#  PRECASTObj <- PRECAST(PRECASTObj, K= 7)

## ----eval = FALSE-------------------------------------------------------------
#  ## backup the fitting results in resList
#  resList <- PRECASTObj@resList
#  PRECASTObj <- selectModel(PRECASTObj)
#  ari_precast <- mclust::adjustedRandIndex(PRECASTObj@resList$cluster[[1]], PRECASTObj@seulist[[1]]$layer_guess_reordered)

## ----eval = FALSE-------------------------------------------------------------
#  seuInt <- PRECASTObj@seulist[[1]]
#  seuInt@meta.data$cluster <- factor(unlist(PRECASTObj@resList$cluster))
#  seuInt@meta.data$batch <- 1
#  seuInt <- Add_embed(PRECASTObj@resList$hZ[[1]], seuInt, embed_name = 'PRECAST')
#  posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
#  seuInt <- Add_embed(posList[[1]], seuInt, embed_name = 'position')
#  Idents(seuInt) <- factor(seuInt@meta.data$cluster)
#  
#  seuInt
#  ## The low-dimensional embeddings obtained by PRECAST are saved in PRECAST reduction slot.

## ----eval = FALSE-------------------------------------------------------------
#  p_sp1 <- SpaPlot(seuInt, item='cluster', point_size = 3, combine = F)[[1]] + cowplot::theme_cowplot() +
#    ggplot2::ggtitle(paste0("PRECAST: ARI=", round(ari_precast, 2)) ) +
#    ggplot2::xlab("row") + ggplot2::ylab("col")
#  seuInt <- AddTSNE(seuInt,n_comp = 2)
#  p_tsne <- dimPlot(seuInt, item='cluster')
#  p_tsne <- p_tsne + cowplot::theme_cowplot() + ggplot2::ggtitle("PRECAST")

## ----eval = FALSE-------------------------------------------------------------
#  seu_drsc <- DR.SC::DR.SC(PRECASTObj@seulist[[1]], K=7, verbose=T)
#  ari_drsc <- mclust::adjustedRandIndex(seu_drsc$spatial.drsc.cluster, PRECASTObj@seulist[[1]]$layer_guess_reordered)
#  p_tsne_drsc <- DR.SC::drscPlot(seu_drsc)
#  p_tsne_drsc <- p_tsne_drsc + ggplot2::ggtitle("DR-SC")
#  p_sp2 <- DR.SC::spatialPlotClusters(seu_drsc)+ cowplot::theme_cowplot()  +
#    ggplot2::ggtitle(paste0("DR-SC ARI=", round(ari_drsc, 2)) )

## ----eval = FALSE, fig.width=10, fig.height=4---------------------------------
#  drawFigs(list(p_sp1, p_sp2), layout.dim = c(1,2))
#  

## ----eval = FALSE, fig.width=10, fig.height=4---------------------------------
#  drawFigs(list(p_tsne, p_tsne_drsc), layout.dim = c(1,2))

## ----eval = FALSE-------------------------------------------------------------
#  sessionInfo()

